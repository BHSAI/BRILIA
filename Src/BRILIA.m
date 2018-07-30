%                              BRILIA
%(B-cell Repertoire Inductive Lineage and Immunosequence Annotator)
%  Written by Donald Lee (dlee@bhsai.org)
%  Developed at BHSAI (www.bhsai.org)
%  Instructions at https://github.com/BHSAI/BRILIA
%
%BASIC USAGE
%  1) Start BRILIA.exe (Win) or BRILIA.sh (Linux) in the command prompt
%
%  2) Annotate input files(s):
%     FMT InputFile(s)                    Param   Value Param  Value Param Value 
%         ------------------------------- ------- ----- ------ ----- ----- ----- 
%       > "C:\My Folder\Grp*_File*.fasta" Species Mouse Strain C57BL Chain H 
%     NOTE: '*' = wildcard string, and '**' = wildcard folder+subfolders
%           See below for Param-Value inputs for BRILIA
%
%  3) Analyze output file(s):
%     runAnalysis: run basic repertoire statistics calculations, such as SHM freq, VDJ usage 
%     plotTree: plot lineage trees for all or certain clusters
%     cmprRep: compare repertoires from multiple files
%
%     FMT Command          OutputFile(s)                     Param Value Param Value
%         -----------      --------------------------------- ----- ----- ----- -----
%       > plotTree         "C:\My Dir\**\Grp*_SeqFile*.csv"
%       > runAnalysis      "C:\My Dir\**\Grp*_SeqFile*.csv"
%       > runGroupAnalysis "C:\My Dir\**\Grp*_SeqFile*.csv"  G1    Grp1  G2    Grp2 
%     NOTE: Type "help [Command]" to view details for each Command's Param-Value inputs.
%
%  *GUIs for BRILIA and analysis tools are:
%       > GUI_BRILIA             %for BRILIA
%       > GUI_plotTree           %for plotting lineage trees
%       > GUI_runAnalysis        %for individual repertoire analysis
%       > GUI_runGroupAnalysis   %for multiple repertoire comparisons
%     
%  PARAM-VALUE INPUTS FOR BRILIA
%
%     Param       Value (defaults = *)     Details
%     ----------- ------------------------ --------------------------------
%     Input       * ""                     Ask user to select input file
%                   "Folder/File*.fa*"     Input sequence files (.fasta/q or .csv). Wildcard string '*' and folder '**'.
%                                          NOTE: The 1st BRILIA input is assumed to be this, so you can omit "Input".
%     Chain       * ""                     Ask user to select IgG chain
%                   H                      Heavy chain
%                   L                      Light chain
%                   HL                     Heavy and Light chains
%     Species     * ""                     Ask user to select database to use from IMGT
%                   human                  human
%                   mouse                  all mouse strain 
%                   macaque                crab-eating macaque
%                   zebrafish              zebrafish
%     Strain      * all                    For mouse only. Use all strains. 
%                   c57bl                  C57BL strains, include C57BL6, C57BL/J
%                   balb                   BALB strains, such as BALB/C
%                   ExactName              NOTE: Some databases are incomplete for certain strains
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
%     SeqRange    * [1,Inf]                Process all sequences 
%                   #                      Process only the #th sequence
%                   [M,N]                  Process Mth to Nth seqeunce (include brackets "[]" , "," , and NO SPACE)
%     CheckSeqDir * y                      Check both fowrad and rev-comp alignment for best annotation
%                   n                      Skip rev-comp alignment (faster if data is pre-processed to + sense only)
%     SkipLineage * n                      Do not skip lineage-based annotation clustering & correction
%                   y                      Skip lineage-based correction, if all sequences are clonally unrelated
%     Output      * ""                     Default output directory at "InputPath/InputName/"
%                   "Out Folder/"          Custom output directory
%                                          Note: If the input is "InputName.fasta", then output is "InputName.BRILIAvN.csv"
%     AutoExit    * y                      Will exit BRILIA local environment when job completes
%                   n                      Will not exit BRILIA local environment

function varargout = BRILIA(varargin)
Version = '4.0.0';
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

SpeciesList = getGeneDatabase('getlist');
ChainList = {'H', 'L', 'HL', 'LH'};
SubFuncNames = [];

P = inputParser;
P.PartialMatching = 1; %Needed for some backward compatibility "Input" and "InputFile" will both work. 
addParameter(P, 'InputFile',     '',      @(x) ischar(x) || iscell(x) || isempty(x));
addParameter(P, 'Chain',         'h',     @(x) ismember({upper(x)}, ChainList));
addParameter(P, 'Species',       '',      @(x) ischar(x) && any(contains(SpeciesList, x, 'ignorecase', true)));
addParameter(P, 'Strain',        'all',   @ischar);
addParameter(P, 'Dgene',         'all',   @(x) ischar(x) && ismember(lower(x), {'all', 'fwd', 'rev'}));
addParameter(P, 'Vgene',         'f',     @(x) ischar(x) && all(ismember(strsplit(lower(x), ','), {'all', 'f', 'p', 'orf'})));
addParameter(P, 'Cores',         'max',   @(x) ischar(x) || isnumeric(x));
addParameter(P, 'CheckSeqDir',   'y',     @(x) ischar(x) && ismember(lower(x), {'y', 'n'}));
addParameter(P, 'StatusHandle',  [],      @(x) ishandle(x) || isempty(x) || strcmpi(class(x), 'matlab.ui.control.UIControl'));
addParameter(P, 'SeqRange',      [1,Inf], @(x) isnumeric(x) || ischar(x));
addParameter(P, 'Resume',        'y',     @(x) ischar(x) && ismember(lower(x), {'y', 'n', 'ask'})); %Resume from an incomplete job
addParameter(P, 'OutputDir',     [],      @(x) ischar(x) || iscell(x) || isempty(x));
addParameter(P, 'BatchSize',     30000,   @(x) isnumeric(x) && x >= 1);
addParameter(P, 'SkipLineage',   'n',     @(x) ischar(x) && ismember(lower(x), {'y', 'n'}));
addParameter(P, 'AutoExit',      'y',     @(x) ischar(x) && ismember(lower(x), {'y', 'n'}));

%Determine the BRILIA run mode
if nargin == 0
    RunInLocalEnv = 1;  %Run in BRILIA local environment.
else
    RunInLocalEnv = 0;  %Run without loading local env. Input provided.
end

while true
    if RunInLocalEnv
        if ~HasShownCredit %Show credits early, once
            showCredits('bhsai', 'imgt');
            HasShownCredit = true;
            fprintf('\nType input commands, ''exit'', or ''help''.  (Wrap files with spaces in quotes, like "C:\\Temp Dir\\)":\n'); %Only show this once
        end
        
        Input = input('BRILIA> ', 's'); %Start with user input
        if strcmpi(Input, 'exit'); return; end
        varargin = cleanCommandLineInput(Input);
        if isempty(varargin); continue; end
    else
        varargin = cleanCommandLineInput(varargin{:});
    end

    %Special parsing of first input
    if ~isempty(varargin) && ischar(varargin{1})
        %Getting version
        if any(strcmpi(varargin{1}, {'ver', 'version', 'getversion'}))
            if RunInLocalEnv
                fprintf('BRILIA VERSION %s\n\n', Version);
                continue
            else
                varargout{1} = Version;
                return
            end
        end
        
        %Getting help
        if any(strcmpi(varargin{1}, {'h', 'i', 'info', 'help', 'showhelp', 'showinfo'}))
            if length(varargin) > 1
                showHelp(varargin{2});
            else
                showHelp('BRILIA');
            end
            if RunInLocalEnv; continue; else; return; end
        end

        %Summon subfunction
        if isempty(SubFuncNames)
            BasicFuncNames = {'system'; 'ls'; 'cd'; 'dir'};
            SubFuncNames = [BasicFuncNames; arrayfun(@(x) x.name(1:end-2), dir(fullfile(findRoot, '**', '*.m')), 'un', 0)];
        end
        if ismember(varargin{1}, SubFuncNames)
            try
                FH = str2func(varargin{1});
                FH(varargin{2:end});
            catch ME
                disp(ME);
            end
            if RunInLocalEnv; continue; else; return; end
        end

        %Check 1st input
        if ~isempty(dir(varargin{1})) || isempty(varargin{1})
            varargin = ['InputFile' varargin]; %#ok<AGROW>
        end
    end

    try
        parse(P, varargin{:});
        P = P.Results;
    catch ME
        fprintf('%s: Error parsing input at line = %d.\n  %s\n', mfilename, ME.stack(1).line, ME.message);
        if RunInLocalEnv; continue; else; return; end
    end

    InputFile = P.InputFile;
    OutputDir = P.OutputDir;
    Chain = strrep(upper(P.Chain), 'LH', 'HL');
    Species = P.Species;
    Strain = P.Strain;
    Dgene = P.Dgene;
    Vgene = P.Vgene;
    Cores = P.Cores;
    BatchSize = round(P.BatchSize);
    SeqRange = round(P.SeqRange);
    Resume = P.Resume;
    StatusHandle = P.StatusHandle;
    CheckSeqDir = P.CheckSeqDir;
    SkipLineage = P.SkipLineage;
    AutoExit = P.AutoExit;

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
            
    OutputFile = cell(size(InputFile));
    for f = 1:length(OutputFile)
        [FilePath, ~, ~, FilePre] = parseFileName(InputFile{f});
        if ~isempty(OutputDir)
            OutputFile{f} = fullfile(OutputDir, FilePre, [FilePre '.BRILIAv' Version(1) '.csv']);
        else
            OutputFile{f} = fullfile( FilePath, FilePre, [FilePre '.BRILIAv' Version(1) '.csv']);
        end
    end

    %==========================================================================
    %BRILIA processing begins 

    DB = getGeneDatabase(Species, 'Strain', Strain, 'Dgene', Dgene, 'Vgene', Vgene);

    showStatus('Setting up parallel computing ...', StatusHandle);
    [~, NumWorkers] = setCores(Cores);
    showStatus(sprintf('  Using %d cores.', NumWorkers), StatusHandle);

    for f = 1:length(InputFile)
        tic
        %----------------------------------------------------------------------
        %File Management

        [OutPath, OutFile, ~, OutFilePre] = parseFileName(OutputFile{f});
        ErrFileName = fullfile(OutPath, [OutFilePre '.Err.csv']);
        TmpFileName = fullfile(OutPath, [OutFilePre '.Tmp.csv']);

        if isdir(OutPath)
            [Success, Msg] = mkdir(OutPath);
            assert(Success, '%s: Could not make the output dir "%s".\n  %s\n', mfilename, OutPath, Msg);
        end
        
        %If there's an incomplete job, Resume from the highest sequence number
        PastTmpFile = arrayfun(@(x) fullfile(x.folder, x.name), dir(fullfile(OutPath, [OutFilePre '*.Tmp.csv'])), 'un', 0);
        PastErrFile = arrayfun(@(x) fullfile(x.folder, x.name), dir(fullfile(OutPath, [OutFilePre  '.Err.csv'])), 'un', 0);
        if ~isempty(PastTmpFile)
            if strcmpi(Resume(1), 'a')
                Choice = lower(input('Found incomplete annotation files. Resume? [y or n] (Enter = y): ', 's'));
                if isempty(Choice) || ~any(strcmpi(Choice(1), {'y', 'n'}))
                    Resume = 'y';
                else
                    Resume = Choice(1);
                end
            end
            
            if strcmpi(Resume(1), 'y')
                showStatus('Resuming from a past incomplete job.', StatusHandle);
                MaxSeqNum = 0;
                if ~isempty(PastErrFile)
                    [VDJdata, ~, ~, ~, Map] = openSeqData(PastErrFile{1});
                    MaxSeqNum = max(cell2mat(VDJdata(:, Map.SeqNum)));
                end
                for b = 1:length(PastTmpFile)
                    [VDJdata, ~, ~, ~, Map] = openSeqData(PastTmpFile{b});
                    MaxSeqNum = max([max(cell2mat(VDJdata(:, Map.SeqNum))), MaxSeqNum]);
                end
                if MaxSeqNum == 0
                    showStatus('No sequences in temp files.', StatusHandle);
                    Resume = 'n'; 
                else
                    SeqRange(1) = MaxSeqNum + 1;
                end
            end
            
            if strcmpi(Resume(1), 'n')
                showStatus('Deleting incomplete job and starting over.\n', StatusHandle);
                delete(PastTmpFile{:});
                delete(PastErrFile{:});
            end
        end
        
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
            [VDJdata, VDJheader, ~, ~, Map] = convertInput2VDJdata(InputFile{f}, 'Chain', Chain, 'SeqRange', SeqRangeB);

            showStatus('Fixing input sequences', StatusHandle);
            [VDJdata, BadLoc] = fixInputSeq(VDJdata, Map);
            KeepLoc(BadLoc) = 0;
            
            %Find potential CDR3 start and end locations using V and J gene seed
            %alignment. Do this here, and not when doing VDJ alignment, because users
            %might have complement sequences which must be flipped.
%           showStatus('Determining sequence direction and CDR3 areas ...', StatusHandle)
%                 VDJdata(KeepLoc, :) = seedCDR3position(VDJdata(KeepLoc, :), Map, DB, 'V', 80,  2, CheckSeqDir);
%                 VDJdata(KeepLoc, :) = seedCDR3position(VDJdata(KeepLoc, :), Map, DB, 'J',  3, 14, 'n');
%                 VDJdata(KeepLoc, :) = seedCDR3position(VDJdata(KeepLoc, :), Map, DB, 'Vk,Vl', 80,  2, CheckSeqDir);
%                 VDJdata(KeepLoc, :) = seedCDR3position(VDJdata(KeepLoc, :), Map, DB, 'Jk,Jl',  3, 14, 'n');

            showStatus('Finding heavy chain VDJ annotations ...', StatusHandle)
            [VDJdata(KeepLoc, :), BadLoc1] = findVDJmatch(VDJdata(KeepLoc, :), Map, DB, 'Update', 'Y');
            
            showStatus('Finding light chain VJ annotations ...', StatusHandle)
            [VDJdata(KeepLoc, :), BadLoc2] = findVJmatch(VDJdata(KeepLoc, :), Map, DB, 'Update', 'Y');
            Idx = find(KeepLoc);
            KeepLoc(Idx(BadLoc1|BadLoc2)) = 0;
            
            showStatus('Fixing indels in V genes ...', StatusHandle);
            VDJdata(KeepLoc, :) = fixGeneIndel(VDJdata(KeepLoc, :), Map, DB);

            showStatus('Accepting F genes over ORF/P ...', StatusHandle);
            VDJdata(KeepLoc, :) = fixDegenVDJ(VDJdata(KeepLoc, :), Map, DB);

            showStatus('Anchoring 104C and 118W/F ...', StatusHandle);
            VDJdata(KeepLoc, :) = constrainGeneVJ(VDJdata(KeepLoc, :), Map, DB);

            showStatus('Moving Nonprod/Invalid/Incomplete Seq to Err file ...', StatusHandle)
            VDJdata = labelSeqQuality(VDJdata, Map, 0.4);
            
            FunctIdx = [Map.hFunct Map.lFunct];
            FunctIdx = FuncIdx(FunctIdx > 0);
            KeepLoc = all(strcmpi(VDJdata(:, FunctIdx), 'Y'), 2);
            if ~all(KeepLoc)
                saveSeqData(ErrFileName, VDJdata(~KeepLoc, :), VDJheader, 'append');
                VDJdata = VDJdata(KeepLoc, :);
                if isempty(VDJdata)
                    warning('%s: No sequences left to annotate in this batch.', mfilename)
                    continue
                end 
            end
            
            if strcmpi(SkipLineage, 'y')
                VDJdata = findCDR(VDJdata, Map, DB, 1:3, 'imgt');
%                 VDJdata = buildVDJalignment(VDJdata, Map, DB);
            end

            %Save remaining sequences to temp dir
            if SeqCount > BatchSize
                CDR3Idx = [Map.hCDR3(1) Map.lCDR3(1)];
                CDR3Idx = CDR3Idx(CDR3Idx > 0);
                CDR3Len = cellfun('length', VDJdata(:, CDR3Idx));

                %Save a separate file per unique CDR3H-L combo
                UnqCDR3Len = unique(CDR3Len, 'rows');
                for k = 1:size(UnqCDR3Len, 1)
                    Idx = ones(size(CDR3Len, 1), 1, 'logical');
                    for q = 1:length(CDR3Idx)
                        Idx = Idx & (CDR3Len(:, q) == UnqCDR3Len(k, q));
                    end
                    if strcmpi(Chain, 'HL')
                        HCDR3Len = UnqCDR3Len(k, 1);
                        LCDR3Len = UnqCDR3Len(k, 2);
                    elseif strcmpi(Chain, 'H')
                        HCDR3Len = UnqCDR3Len(k, 1);
                        LCDR3Len = 0;
                    elseif strcmpi(Chain, 'L')
                        HCDR3Len = 0;
                        LCDR3Len = UnqCDR3Len(k, 1);
                    end
                    SaveName = [TempDir 'CDR3-' num2str(HCDR3Len) '-' num2str(LCDR3Len) '.Raw.csv'];
                    saveSeqData(SaveName, VDJdata(Idx, :), VDJheader, 'append');
                end
            else %If file is small, just save in one file.    
                SaveName = [TempDir 'Raw.csv'];
                saveSeqData(SaveName, VDJdata, VDJheader, 'append');
            end
        end

        %Part 2 does everything AFTER intial annotations
        GrpNumAdj = 0;
        FileList = dir(fullfile(TempDir, '*Raw.csv'));

        if strcmpi(SkipLineage, 'n')
            for t = 1:length(FileList)
                %Check if a Final.csv file already exists to skip
                TempRawFileName = FileList(t).name;
                TempOutFileName = strrep(TempRawFileName, 'Raw.csv', 'Final.csv');
                if exist(fullfile(TempDir, TempOutFileName), 'file'); continue; end
                
                %Reload the sequence file with same-length CDR3s
                showStatus(sprintf('Processing %s ...', TempRawFileName), StatusHandle);
                [VDJdata, VDJheader, ~, ~, Map] = openSeqData(fullfile(TempDir, TempRawFileName));
                if size(VDJdata, 1) < 1; continue; end
                
                %Cluster the data based variable region and hamming dist of DevPerc%.
                showStatus('Clustering by lineage ...', StatusHandle)
                VDJdata = clusterByJunction(VDJdata, Map);
                VDJdata = spliceData(VDJdata, Map); %Want to make it parfor-capable
                VDJdata = clusterByLineage(VDJdata, Map, 'shmham');

                %Set all groups to have same annotation and VMDNJ lengths.
                showStatus('Correcting annotations by lineage ...', StatusHandle)
                VDJdata = conformGeneGroup(VDJdata, Map, DB);

                %Get better D match based on location of consensus V J mismatches.
                showStatus('Refining D annotations ...', StatusHandle)
                VDJdata = findBetterD(VDJdata, Map, DB);

                %Trim V, D, J edges and extract better N regions
                showStatus('Trimming N regions ...', StatusHandle)
                VDJdata = trimGeneEdge(VDJdata, Map, DB);
                
                VDJdata = joinData(VDJdata, Map);

                %Finalize VDJdata details and CDR 1, 2, 3 info
                VDJdata = padtrimSeqGroup(VDJdata, Map, 'grpnum', 'trim', 'Seq'); %will only remove "x" before and after Seq if they all have it. 
                VDJdata = findCDR1(VDJdata, Map, DB);
                VDJdata = findCDR2(VDJdata, Map, DB);
                VDJdata = findCDR3(VDJdata, Map, DB, 'IMGT'); %removes the 104C and 118W from CDR3, and adjusts the CDR3 length to true IMGT length
                VDJdata = buildVDJalignment(VDJdata, Map, DB);
                
                %Renumbering groups since they're done in batches
                GrpNums = cell2mat(VDJdata(:, Map.GrpNum)) + GrpNumAdj;
                GrpNumAdj = max(GrpNums); 
                VDJdata(:, Map.GrpNum) = num2cell(GrpNums);

                saveSeqData([TempDir TempOutFileName], VDJdata, VDJheader, 'append');
            end
            %======================================================================
            %Move folder to the correct destination folder

            %Combine final annotations to a single file, then move to destination
            FinalFileList = arrayfun(@(x) fullfile(TempDir, x.name), dir(fullfile(TempDir, '*Final.csv')), 'unif', false);
            PrevFinalFile = fullfile(TempDir, OutFile);
            if exist(PrevFinalFile, 'file') %If a real output file was interrupted and left in Temp folder, delete.
                try
                    delete(PrevFinalFile);
                catch
                    warning('%s: Could not delete incomplete final file "%s".', mfilename, PrevFinalFile);
                end
            end
            if length(FinalFileList) > 1
                SrcFile = fullfile(TempDir, OutFile);
                DstFile = fullfile(OutPath, OutFile);
                combineSeqData(FinalFileList, SrcFile);
                try
                    movefile(SrcFile, DstFile, 'f');
                    delete(FinalFileList{:});
                catch
                    warning('%s: Could not move Final files from "%s" to destination "%s".', mfilename, SrcFile, DstFile);
                end
            elseif length(FinalFileList) == 1
                SrcFile = FinalFileList{1};
                DstFile = fullfile(OutPath, OutFile);
                try
                    movefile(SrcFile, DstFile, 'f');
                catch
                    warning('%s: Could not move Final files from "%s" to destination "%s".', mfilename, SrcFile, DstFile);
                end
            end
        end

        %Combine raw annotations to a single file, then move to destination
        RawFileList = arrayfun(@(x) fullfile(TempDir, x.name), dir(fullfile(TempDir, '*Raw.csv')), 'unif', false);
        if length(RawFileList) > 1
            [~, RawFileNameOnly] = parseFileName(TmpFileName);
            SrcFile = TmpFileName;
            DstFile = fullfile(OutPath, RawFileNameOnly);
            try
                combineSeqData(RawFileList, SrcFile);
                movefile(SrcFile, DstFile);
                delete(RawFileList{:});
            catch
                warning('%s: Could not move Raw files from "%s" to destination "%s".', mfilename, TempDir, DstFile);
            end
        elseif length(RawFileList) == 1
            [~, RawFileNameOnly] = parseFileName(TmpFileName);
            SrcFile = RawFileList{1};
            DstFile = fullfile(OutPath, RawFileNameOnly);
            try
                movefile(SrcFile, DstFile, 'f');
            catch
                warning('%s: Could not move Raw files from "%s" to destination "%s".', mfilename, SrcFile, DstFile);
            end
        end

        %Move the err file out
        if exist(ErrFileName, 'file')
            [~, ErrFileNameOnly] = parseFileName(ErrFileName);
            SrcFile = ErrFileName;
            DstFile = fullfile(OutPath, ErrFileNameOnly);
            try
                movefile(SrcFile, DstFile, 'f');
            catch
                warning('%s: Could not move Err files from "%s" to destination "%s".', mfilename, SrcFile, DstFile);
            end
        end

        prepTempDir(TempDir, 'delete');
        showStatus(sprintf('Finished in %0.1f sec.', toc), StatusHandle);
    end
    varargout{1} = OutputFile;
    
    if RunInLocalEnv && strcmpi(AutoExit, 'n'); continue; else; return; end
end