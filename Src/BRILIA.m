%                              BRILIA
%(B-cell Repertoire Inductive Lineage and Immunosequence Annotator)
%  Written by Donald Lee (dlee@bhsai.org)
%  Developed at BHSAI (www.bhsai.org)
%  Instructions at https://github.com/BHSAI/BRILIA
%
%BASIC USAGE
%  1) Start BRILIA.exe (Win) or run_BRILIA.sh (Linux) in the command prompt
%
%  2) Annotate input file(s):
%        InputFile(s)                    Param   Value Param  Value Param Value 
%        ------------------------------- ------- ----- ------ ----- ----- -----
%     >  "C:\My Folder\Grp*_File*.fasta" Species Mouse Strain C57BL Chain H 
%     
%     NOTE: '*'  = wildcard string
%           '**' = wildcard folder + subfolders
%           See below for Param-Value inputs for BRILIA
%
%  3) Analyze output file(s):
%     plotTree: plot lineage trees for all or certain clusters
%     runAnalysis: run basic repertoire statistics calculations, such as SHM freq, VDJ usage 
%     runGroupAnalysis: compares repertoires from multiple files
%
%         Command          OutputFile(s)                     Param Value Param Value
%         -----------      --------------------------------- ----- ----- ----- -----
%       > plotTree         "C:\My Dir\**\Grp*_SeqFile*.csv"
%       > runAnalysis      "C:\My Dir\**\Grp*_SeqFile*.csv"
%       > runGroupAnalysis "C:\My Dir\**\Grp*_SeqFile*.csv"  G1    Grp1  G2    Grp2 
%
%     NOTE: Type "help [Command]" to view details
%
%  4) Analyze output file(s) via GUI:
%       > GUI_BRILIA             %for BRILIA
%       > GUI_plotTree           %for plotting lineage trees
%       > GUI_runAnalysis        %for individual repertoire analysis
%       > GUI_runGroupAnalysis   %for multiple repertoire comparisons
%     
%  PARAM-VALUE INPUTS FOR BRILIA
%
%     Param       Value (* = default)      Details
%     ----------- ------------------------ --------------------------------
%     Input       * []                     Ask user to select input file
%                   "Folder\File*.fa*"     Input sequence files (.fa* or .csv). * = any string, ** = any subfolder.
%                                          NOTE: The 1st BRILIA input is assumed to be this, so you can omit "Input".
%     OutputDir   * []                     Default output directory at "InputPath\InputName\"
%                   "OutDir\"              Custom output directory. Note: Cannot specify OutputFile AND OutputDir.
%     OuptutFile  * []                     Default output file "InputPath\InputName\InputName.BRILIAvN.csv"
%                   "OutDir\OutFile.csv"   Custom output CSV file. All temp files will placed in same folder
%     Chain       * []                     Ask user to select IgG chain
%                   h                      Heavy chain
%                   l                      Light chain
%                   hl                     Heavy and Light chains
%     Cutoff      * 0                      Do not merge similar sequences
%                   F   (0 to 0.99)        Merge similar sequences that are fraction 0 < F < 1 of sequence length
%                   N   (integer >= 1)     Merge similar sequences that are integer N >= 1 Hamming distance
%     Species     * []                     Ask user to select database to use from IMGT
%                   human                  human
%                   mouse                  all mouse strain 
%                   macaque                crab-eating macaque
%                   zebrafish              zebrafish
%                   "Species"              species name in the Database folder where BRILIA exec file is.
%                                          NOTE: you can add custom database by adding a folder with IGXX.fa files
%     Strain      * all                    For mouse only. Use all strains. 
%                   c57bl                  C57BL strains, include C57BL6, C57BL/J
%                   balb                   BALB strains, such as BALB/C
%                   "Strain"               The prefix of strains that match
%     Dgene       * all                    Foward and inverse are okay
%                   fwd                    Foward only
%                   inv                    Inverse only
%     Vgene       * f                      Functional only
%                   p                      Psueudo genes only
%                   orf                    Open reading frame only
%                   all                    All of the above
%     Batch       * 30000                  Process 30000 sequences per batch to prevent memory overload
%                   #                      Process # sequences per batch
%     Cores       * max                    Use maximum number of cores
%                   #                      Use # number of cores
%     Resume      * y                      Resume from an interrupted job 
%                   n                      Do not resume from an interrupted job, ex server outage
%                   ask                    Ask user to confirm to resume if temp files are found
%     SeqRange    * [1,inf]                Process all sequences 
%                   #                      Process only the #th sequence
%                   [M,N]                  Process Mth to Nth seqeunce (include brackets "[]" , "," , and NO SPACE)
%     MinQuality  * 2                      Min Phred Score (ASCII 33). DNA bases w/ read error >= 2% will be "N". 
%                   char                   Min Phred Score (ASCII 33). 
%                                          NOTE: this is only for pre-processing fastq files.
%     CheckSeqDir * y                      Check both fowrad and rev-comp alignment for best annotation
%                   n                      Skip rev-comp alignment (faster if data is pre-processed to + sense only)
%     SkipLineage * n                      Do not skip lineage-based annotation clustering & correction
%                   y                      Skip lineage-based correction, if all sequences are clonally unrelated
%     AutoExit    * n                      Do not exit BRILIA local environment when job completes
%                   y                      Exit BRILIA local environment when job completes

function varargout = BRILIA(varargin)
Version = '3.5.1'; 
varargout = cell(1, nargout);
HasShownCredit = false;

%--------------------------------------------------------------------------
%For running in matlab, make sure BRILIA paths are added correctly
if ~isdeployed
    MainPath = fileparts(fileparts(mfilename('fullpath')));
    if ~any(strcmpi(strsplit(path, pathsep), MainPath))
        if strcmpi(input('Add BRILIA path to matlab? y or n: ', 's'), 'y')
            fprintf('Adding BRILIA path to Matlab.\n');
            addpath(genpath(MainPath));
        else
            fprintf('Aborting BRILIA without adding to matlab path.\n');
            return
        end
    end
end

%--------------------------------------------------------------------------
%Handle various inputs from matlab or OS command lines
RunInLocalEnv = (nargin == 0); % 1 (run locally in BRILIA's own session?)
SpeciesList = getGeneDatabase('getlist');
ChainList = {'H', 'L', 'HL', 'LH'};
SubFuncNames = [];

P = inputParser; 
addParameter(P, 'InputFile',     '',      @(x) ischar(x) || iscell(x) || isempty(x));
addParameter(P, 'OutputDir',     '',      @(x) ischar(x) || iscell(x) || isempty(x));
addParameter(P, 'OutputFile',    '',      @(x) ischar(x) || iscell(x) || isempty(x));
addParameter(P, 'Chain',         'h',     @(x) ismember({upper(x)}, ChainList));
addParameter(P, 'Cutoff',        0,       @(x) isnumeric(x) && x >= 0);
addParameter(P, 'Species',       '',      @(x) ischar(x) && any(contains(SpeciesList, x, 'ignorecase', true)));
addParameter(P, 'Strain',        'all',   @ischar);
addParameter(P, 'Dgene',         'all',   @(x) ischar(x) && ismember(lower(x), {'all', 'fwd', 'rev'}));
addParameter(P, 'Vgene',         'f',     @(x) ischar(x) && all(ismember(strsplit(lower(x), ','), {'all', 'f', 'p', 'orf'})));
addParameter(P, 'BatchSize',     30000,   @(x) isnumeric(x) && x >= 1);
addParameter(P, 'Cores',         'max',   @(x) ischar(x) || isnumeric(x));
addParameter(P, 'SeqRange',      [1,Inf], @(x) isnumeric(x) || ischar(x));
addParameter(P, 'MinQuality',    '2',     @(x) ischar(x) || isnumeric(x)); %ASCII_BASE=33, '2' = P_error 0.01995
addParameter(P, 'StatusHandle',  [],      @(x) ishandle(x) || isempty(x) || strcmpi(class(x), 'matlab.ui.control.UIControl'));
addParameter(P, 'Resume',        'y',     @(x) ischar(x) && ismember(lower(x), {'y', 'n', 'ask'})); %Resume from an incomplete job
addParameter(P, 'CheckSeqDir',   'y',     @(x) ischar(x) && ismember(lower(x), {'y', 'n'}));
addParameter(P, 'SkipLineage',   'n',     @(x) ischar(x) && ismember(lower(x), {'y', 'n'}));
addParameter(P, 'AutoExit',      'n',     @(x) ischar(x) && ismember(lower(x), {'y', 'n'}));
addParameter(P, 'SettingFile',   '',      @(x) isempty(x) || ~isempty(dir(x))); %This is kept for backward compatibility only. Will be removed.

while true
    if RunInLocalEnv
        if ~HasShownCredit %Show credits early, once
            showCredits('bhsai', 'imgt');
            HasShownCredit = true;
            fprintf('\nType input commands, ''exit'', or ''help''.  (Wrap files with spaces in quotes, like "C:\\Temp Dir\\)":\n'); %Only show this once
        end
        
        Input = input('BRILIA> ', 's'); %Start with user input
        if isempty(Input) 
            continue
        elseif strcmpi(Input, 'exit')
            return
        end
        varargin = cleanCommandLineInput(Input);
    else
        varargin = cleanCommandLineInput(varargin{:});
    end

    %Special parsing of first input
    if ~isempty(varargin) && ischar(varargin{1})
        if any(strcmpi(varargin{1}, {'ver', 'version', 'getversion'}))
            if RunInLocalEnv
                fprintf('BRILIA VERSION %s\n\n', Version);
                continue
            else
                varargout{1} = Version;
                return
            end
        end
        
        if any(strcmpi(varargin{1}, {'?', 'h', 'i', 'info', 'help', 'showhelp', 'showinfo'}))
            if length(varargin) > 1
                showHelp(varargin{2});
            else
                showHelp('BRILIA');
            end
            if RunInLocalEnv; continue; else; return; end
        end
        
        %if user has inputed "help list", show ALL available functions

        if isempty(SubFuncNames)
            BasicFuncNames = {'system'; 'ls'; 'cd'; 'dir'};
            SubFuncNames = [BasicFuncNames; arrayfun(@(x) x.name(1:end-2), dir(fullfile(findRoot, '**', '*.m')), 'un', 0)];
        end
        if ismember(varargin{1}, SubFuncNames)
            try
                FH = str2func(varargin{1});
                FH(varargin{2:end});
            catch ME
                disp(ME)
                disp(varargin)
                for j = 1:length(ME.stack)
                    disp(ME.stack(j))
                end
            end
            if RunInLocalEnv; continue; else; return; end
        end

        if ~isempty(dir(varargin{1})) || isempty(varargin{1})
            varargin = ['InputFile' varargin]; %#ok<AGROW>
        end
    end

    try
        %Backward compatibility - change Vfunction to Vgene, Ddirection to Dgene, and ignore DevPerc
        CharLoc = cellfun('isclass', varargin, 'char');
        varargin(CharLoc) = regexprep(varargin(CharLoc), {'Vfunction', 'Ddirection'}, {'Vgene', 'Dgene'}, 'ignorecase');
        
        [Ps, ~, ReturnThis] = parseInput(P, varargin{:}); %Default from parseInput: Partial Matching, Case-Insensitive, Keep Unmatched, Struct Expand
        if ReturnThis && ~RunInLocalEnv
            varargout{1} = Ps;
            return
        end
        
        %Override defaults with what is in the SettingFile
        if ~isempty(Ps.SettingFile)
            Ps = readSettingFile(Ps.SettingFile, Ps);
        end

    catch ME
        if RunInLocalEnv 
            fprintf('%s: Error parsing input at line = %d.\n  %s\n', mfilename, ME.stack(1).line, ME.message);
            continue
        else
            error('%s: Error parsing input at line = %d.\n  %s\n', mfilename, ME.stack(1).line, ME.message);
        end
    end
    
    InputFile = Ps.InputFile;
    OutputDir = Ps.OutputDir;
    OutputFile = Ps.OutputFile;
    Chain = strrep(upper(Ps.Chain), 'LH', 'HL');
    Cutoff = Ps.Cutoff;
    Species = Ps.Species;
    Strain = Ps.Strain;
    Vgene = Ps.Vgene;
    Dgene = Ps.Dgene;
    Cores = Ps.Cores;
    BatchSize = round(Ps.BatchSize);
    SeqRangeT = round(Ps.SeqRange);
    StatusHandle = Ps.StatusHandle;
    Resume = Ps.Resume;
    CheckSeqDir = Ps.CheckSeqDir;
    SkipLineage = Ps.SkipLineage;
    AutoExit = Ps.AutoExit;
    MinQuality = Ps.MinQuality;

    if ~HasShownCredit
        showCredits('bhsai', 'imgt');
        HasShownCredit = true;
    end
    
    if isempty(Chain)
        fprintf('Error: Did not specify the "Chain" parameter. Valid values are:\n');
        fprintf('  %s\n', ChainList{:});
        if RunInLocalEnv; continue; else; return; end
    end
    
    if isempty(Species)
        fprintf('Error: Did not specify the "Species" parameter. Valid values are:\n');
        fprintf('  %s\n', SpeciesList{:});
        if RunInLocalEnv; continue; else; return; end
    end
    
    %--------------------------------------------------------------------------
    %Input and Output Files

    if isempty(InputFile) %Ask user to choose
        InputFile = openFileDialog('*.fa*;*.*sv', 'Select the input sequence files', 'multiselect', 'on');
    elseif ischar(InputFile) %Search for all files that matches
        InputFile = dir(InputFile);
        InputFile = fullfile({InputFile.folder}, {InputFile.name});
    end
    
    if isempty(InputFile)
        fprintf('No valid file was selected.\n\n');
        if RunInLocalEnv; continue; else; return; end
    end

    %Correct chain from HL to H if dealing with fasta/q
    if strcmpi(Chain, 'HL') && any(endsWith(InputFile, {'.fa', '.fasta', '.fastq'}, 'ignorecase', true))
        fprintf('For Chain = HL option, only delimited (ie, *.csv) files can be used.\n')
        fprintf('The delimited file must have defined "H-Seq" and "L-Seq" columns.\n'); 
        if RunInLocalEnv; continue; else; return; end
    end
    
    
    %Make sure # of output files = # of input files
    if ~isempty(OutputFile)
        if ~isempty(OutputDir)
            if RunInLocalEnv
                fprintf('%s: Cannot use inputs for OutputDir AND OuputFile. Pick one.\n', mfilename);
                continue
            else
                error('%s: Cannot use inputs for OutputDir AND OuputFile. Pick one.', mfilename);
            end            
        end
        if ischar(OutputFile)
            OutputFile = {OutputFile};
        end
        if length(OutputFile) ~= length(InputFile)
            if RunInLocalEnv
                fprintf('%s: Mismatched number of input (%d) and output files (%d).\n', mfilename, length(InputFile), length(OutputFile));
                continue
            else
                error('%s: Mismatched number of input (%d) and output files (%d).', mfilename, length(InputFile), length(OutputFile));
            end            
        end
        if ~RunInLocalEnv && nargout >= 1
            varargout{1} = OutputFile;
        end
    else
        OutputFile = cell(size(InputFile));
        for f = 1:length(OutputFile)
            [FilePath, ~, ~, FilePre] = parseFileName(InputFile{f});
            if ~isempty(OutputDir)
                OutputFile{f} = fullfile(OutputDir, FilePre, [FilePre '.BRILIAv' Version(1) '.csv']);
            else
                OutputFile{f} = fullfile(FilePath,  FilePre, [FilePre '.BRILIAv' Version(1) '.csv']);
            end
            if ~RunInLocalEnv && nargout >= 1
                varargout{1} = OutputFile;
            end
        end
    end

    %==========================================================================
    %BRILIA processing begins 

    DB = getGeneDatabase(Species, 'Strain', Strain, 'Dgene', Dgene, 'Vgene', Vgene);

    showStatus('Setting up parallel computing ...', StatusHandle);
    [~, NumWorkers] = setCores(Cores);
    showStatus(sprintf('  Using %d cores.', NumWorkers), StatusHandle);

    for f = 1:length(InputFile)
        SeqRange = SeqRangeT; %Reset this every time, otherwise you'll get odd results.
        
        TicInputFile = tic;
        %------------------------------------------------------------------
        %File Management

        [OutPath, ~, ~, OutFilePre] = parseFileName(OutputFile{f});
        if ~isdir(OutPath)
            [Success, Msg] = mkdir(OutPath);
            if ~Success 
                warning('%s: Could not make the output dir "%s".\n  %s\n  Skipping file.\', mfilename, OutPath, Msg);
                continue
            end
        end
        ErrFileName = fullfile(OutPath, [OutFilePre '.Err.csv']);
        TmpFileName = fullfile(OutPath, [OutFilePre '.Tmp.csv']);
        RawFileName = fullfile(OutPath, [OutFilePre '.Raw.csv']);
        TmpExist = ~isempty(dir(TmpFileName));
        ErrExist = ~isempty(dir(ErrFileName));
        RawExist = ~isempty(dir(ErrFileName));
        
        %Resume from the highest sequence number
        if strcmpi(Resume, 'n')
            showStatus('Deleting incomplete job and starting over.\n', StatusHandle);
            if TmpExist; delete(TmpFileName); end
            if ErrExist; delete(ErrFileName); end
            if RawExist; delete(RawFileName); end
            TmpExist = 0;
            ErrExist = 0;
            RawExist = 0;
        end
        
        %Resume from the highest sequence number
        if TmpExist
            if strcmpi(Resume, 'a')
                Choice = input('Found incomplete annotation files. Resume? [y or n] (Enter = y): ', 's');
                if isempty(Choice) || ~ismember(lower(Choice(1)), {'y', 'n'})
                    Resume = 'y';
                else
                    Resume = Choice(1);
                end
            end

            if strcmpi(Resume, 'y')
                showStatus('Resuming from a past incomplete job.', StatusHandle);
                MaxSeqNum = 0;
                if ErrExist
                    [VDJdata, ~, ~, ~, Map] = openSeqData(ErrFileName);
                    MaxSeqNum = max(cell2mat(VDJdata(:, Map.SeqNum)));
                end
                [VDJdata, ~, ~, ~, Map] = openSeqData(TmpFileName);
                MaxSeqNum = max([max(cell2mat(VDJdata(:, Map.SeqNum))), MaxSeqNum]);
                if MaxSeqNum == 0
                    showStatus('No sequences in temp files.', StatusHandle);
                    Resume = 'n'; 
                else
                    SeqRange(1) = MaxSeqNum + 1;
                end
            end
        end
        
        if ~RawExist
            %Set the sequence range
            if numel(SeqRange) == 1
                SeqRange = repelem(SeqRange, 1, 2);
            end
            MaxSeqCount = countSeq(InputFile{f});
            SeqRange(SeqRange > MaxSeqCount) = MaxSeqCount;
            SeqRange(SeqRange < 1) = 1;
            SeqCount = diff(SeqRange) + 1;

            %Part 1: performs V(D)J annotation
            showStatus(sprintf('Opening "%s" ...', InputFile{f}), StatusHandle); 
            for b = 1:ceil(SeqCount/BatchSize)
                SeqRangeB = [SeqRange(1)+BatchSize*(b-1)  SeqRange(1)+b*BatchSize-1]; %batch seq range
                SeqRangeB(end) = min([SeqRangeB(end) SeqCount]);
                KeepLoc = ones(diff(SeqRangeB)+1, 1, 'logical');

                showStatus(sprintf('Processing sequences %d to %d (out of %d) ...', SeqRangeB(1), SeqRangeB(2), SeqCount), StatusHandle);
                [VDJdata, VDJheader, ~, ~, Map] = convertInput2VDJdata(InputFile{f}, 'Chain', Chain, 'SeqRange', SeqRangeB, 'MinQuality', MinQuality);

                showStatus('Fixing input sequences', StatusHandle);
                [VDJdata, BadLoc] = fixInputSeq(VDJdata, Map);
                KeepLoc(BadLoc) = 0;

                showStatus('Determining sequence direction and CDR3 areas ...', StatusHandle)
                VDJdata(KeepLoc, :) = seedCDR3position(VDJdata(KeepLoc, :), Map, DB, 'V',     40,  2, CheckSeqDir);
                VDJdata(KeepLoc, :) = seedCDR3position(VDJdata(KeepLoc, :), Map, DB, 'J',      3, 14, 'n');
                VDJdata(KeepLoc, :) = seedCDR3position(VDJdata(KeepLoc, :), Map, DB, 'Vk,Vl', 40,  2, CheckSeqDir);
                VDJdata(KeepLoc, :) = seedCDR3position(VDJdata(KeepLoc, :), Map, DB, 'Jk,Jl',  3, 14, 'n');

                showStatus('Finding heavy chain VDJ annotations ...', StatusHandle)
                [VDJdata(KeepLoc, :), BadLoc1] = findVDJmatch(VDJdata(KeepLoc, :), Map, DB, 'Update', 'Y');

                showStatus('Finding light chain VJ annotations ...', StatusHandle)
                [VDJdata(KeepLoc, :), BadLoc2] = findVJmatch(VDJdata(KeepLoc, :), Map, DB, 'Update', 'Y');
                Loc = find(KeepLoc);
                KeepLoc(Loc(BadLoc1|BadLoc2)) = 0;

                showStatus('Fixing insertions/deletions in V genes ...', StatusHandle);
                VDJdata(KeepLoc, :) = fixGeneIndel(VDJdata(KeepLoc, :), Map, DB);

                showStatus('Accepting F genes instead of ORF/P ...', StatusHandle);
                VDJdata(KeepLoc, :) = fixDegenVDJ(VDJdata(KeepLoc, :), Map, DB);

                showStatus('Anchoring 104C and 118W/F ...', StatusHandle);
                VDJdata(KeepLoc, :) = constrainGeneVJ(VDJdata(KeepLoc, :), Map, DB);

                showStatus('Moving non-functional Seq to Err file ...', StatusHandle);
                VDJdata = labelSeqQuality(VDJdata, Map, 0.4);

                FunctIdx = nonzeros([Map.hFunct Map.lFunct]);
                KeepLoc = all(strcmpi(VDJdata(:, FunctIdx), 'Y'), 2);
                if ~all(KeepLoc)
                    saveSeqData(ErrFileName, VDJdata(~KeepLoc, :), VDJheader, 'append');
                    VDJdata = VDJdata(KeepLoc, :);
                    if isempty(VDJdata)
                        showStatus(sprintf('No sequences left to annotate in batch #%d.', b), StatusHandle);
                        continue
                    end 
                end

                if strcmpi(SkipLineage, 'y')
                    VDJdata = findCDR(VDJdata, Map, DB, 1:3, 'imgt');
                    VDJdata = buildVDJalignment(VDJdata, Map, DB);
                end

                saveSeqData(TmpFileName, VDJdata, VDJheader, 'append');
            end
            if exist(TmpFileName, 'file')
                movefile(TmpFileName, RawFileName);
            else
                showStatus(sprintf('There was no functional sequences left.'), StatusHandle);
                continue
            end
        else
            showStatus(sprintf('Have the annotation file already.'), StatusHandle);
        end

        %Part 2 does everything AFTER intial annotations
        if strcmpi(SkipLineage, 'n')
            showStatus(sprintf('Processing %s ...', RawFileName), StatusHandle);
            [VDJdata, VDJheader, ~, ~, Map] = openSeqData(RawFileName);
            if isempty(VDJdata); continue; end
    
%CODING_NOTE: This will be added for future release due to HTS paired end reads with poor quality edge reads.
%Still under development and testing with other datasets.
%2018-09-11
%             showStatus('Cleaning 5'' and 3'' mismatched ends ...', StatusHandle);
%             VDJdata = cleanSeqEnds(VDJdata, Map);
            
            showStatus('Clustering by lineage ...', StatusHandle)
            VDJdata = clusterByJunction(VDJdata, Map);
            VDJdata = spliceData(VDJdata, Map); %Make it parfor-capable
            VDJdata = clusterByLineage(VDJdata, Map, 'shmham');

            showStatus('Correcting annotations by lineage ...', StatusHandle)
            VDJdata = conformGeneGroup(VDJdata, Map, DB);     
            
            showStatus('Refining D annotations ...', StatusHandle)
            VDJdata = findBetterD(VDJdata, Map, DB);

            showStatus('Trimming N regions ...', StatusHandle)
            VDJdata = trimGeneEdge(VDJdata, Map, DB);
                
            showStatus('Removing sequences that are too similar ...', StatusHandle);
            VDJdata = mergeSimilarSeq(VDJdata, Map, Cutoff);

            VDJdata = joinData(VDJdata, Map);

            showStatus('Finalizing annotations ...', StatusHandle);
            VDJdata = padtrimSeqGroup(VDJdata, Map, 'grpnum', 'trim', 'Seq'); 
            VDJdata = findCDR(VDJdata, Map, DB, 1:3, 'imgt');
            VDJdata = buildVDJalignment(VDJdata, Map, DB);
                
            saveSeqData(OutputFile{f}, VDJdata, VDJheader);
        end
        showStatus(sprintf('Finished in %0.1f sec.', toc(TicInputFile)), StatusHandle);
    end
    
    if RunInLocalEnv && strcmpi(AutoExit, 'n'); continue; else; return; end
end