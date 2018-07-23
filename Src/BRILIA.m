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
%     Resume      * n                      Do not resume from an interrupted job, ex server outage
%                   y                      Resume from an interrupted job 
%     Range       * [1,Inf]                Process all sequences 
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
ChainList = {'H', 'L', 'HL'};

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
addParameter(P, 'Range',         [1,Inf], @(x) isnumeric(x) || ischar(x));
addParameter(P, 'Resume',        'n',     @(x) ischar(x) && ismember(lower(x), {'y', 'n'})); %Resumes from known raw file
%addParameter(P, 'ResumeFrom',    '',      @(x) ischar(x)); %File directory storing the *Raw.csv file(s).
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
        if ~HasShownCredit
            showCredits('bhsai', 'imgt');
            HasShownCredit = true;
            fprintf('\nType input commands, ''exit'', or ''help''.  (Wrap files with spaces in quotes, like "C:\\Temp Dir\\)":\n'); %Only show this once
        end
        
        Input = input('BRILIA> ', 's'); %Start with user input
        if strcmpi(Input, 'exit'); return; end
        varargin = cleanCommandLineInput(Input);
        if isempty(varargin); continue; end
    else
        %Will show credits AFTER ensuring BRILIA annotation will run
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
        
        %Getting help for this
        if any(strcmpi(varargin{1}, {'h', 'i', 'info', 'help', 'showhelp', 'showinfo'}))
            if length(varargin) > 1
                showHelp(varargin{2});
            else
                showHelp('BRILIA');
            end
            if RunInLocalEnv; continue; else; return; end
        end

        %Redirect to BRILIA subfunction
        SubFuncNames = {'pwd', 'cd', 'dir', 'ls', ...
                        'plotTree', 'runAnalysis', 'cmprRep', ...
                        'setCores', 'showHelp', ...
                        'GUI_BRILIA', 'GUI_plotTree'}; %For security, only accept allowed function calls!
        if ismember(varargin{1}, SubFuncNames)
            try
                FH = str2func(varargin{1});
                FH(varargin{2:end});
            catch ME
                disp(ME);
            end
            if RunInLocalEnv; continue; else; return; end
        end

        %Check it's a file or correct varargin for parm-value parsing later.
        if exist(varargin{1}, 'file') || isempty(varargin{1})
            varargin = ['InputFile' varargin]; %#ok<AGROW>
        end
    end

    %Parse the inputs and return if this was a code-2-code summon
    try
        parse(P, varargin{:});
        P = P.Results;
    catch ME
        fprintf('%s: Error parsing input at line = %d.\n', mfilename, ME.stack(1).line);
        if RunInLocalEnv; continue; else; return; end
    end

    %Show credit for ~RunInLocalEnv if it hasn't show it yet
    if ~HasShownCredit
        showCredits('bhsai', 'imgt');
        HasShownCredit = true;
    end

    InputFile = P.InputFile;
    UserOutDir = P.OutputDir;
    Chain = upper(P.Chain);
    Species = P.Species;
    Strain = P.Strain;
    Ddirection = P.Dgene;
    Vfunction = P.Vgene;
    Cores = P.Cores;
    CheckSeqDir = P.CheckSeqDir;
    Resume = P.Resume;
    SeqRange = P.Range;
    StatusHandle = P.StatusHandle;
    BatchSize = P.BatchSize;
    SkipLineage = P.SkipLineage;
    AutoExit = P.AutoExit;

    %Make sure Chain and Species are provided. These cannot be empty.
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

    %Get the full file names of input sequences
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

    %Determine output file names
    OutputFile = cell(size(InputFile));
    for f = 1:length(OutputFile)
        [OutPath, ~, ~, OutFilePre] = parseFileName(InputFile{f});
        if ~isempty(UserOutDir)
            OutputFile{f} = fullfile(UserOutDir, OutFilePre, [OutFilePre '.BRILIAv' Version(1) '.csv']);
        else
            OutputFile{f} = fullfile(OutPath, OutFilePre, [OutFilePre '.BRILIAv' Version(1) '.csv']);
        end
    end

    %==========================================================================
    %BRILIA processing begins 

    DB = getGeneDatabase(Species);
    DB = filterGeneDatabase(DB, 'Strain', Strain, 'Ddirection', Ddirection, 'Vfunction', Vfunction);

    showStatus('Setting up parallel computing ...', StatusHandle);
    [~, NumWorkers] = setCores(Cores);
    showStatus(sprintf('  Using %d cores.', NumWorkers), StatusHandle);

    for f = 1:length(InputFile)
        tic
        %----------------------------------------------------------------------
        %File Management

        %Specify the Temp folder and Raw and Err annotation files
        [OutPath, OutFile, ~, OutFilePre] = parseFileName(OutputFile{f}, 'ignorefilecheck');
        TempDir = fullfile(OutPath, 'Temp', filesep);
        ErrFileName = fullfile(TempDir, [OutFilePre '.Err.csv']);
        RawFileName = fullfile(TempDir, [OutFilePre '.Raw.csv']);

%         %If using ResumeFrom, try to move Raw.csv into the TempDir.
%         if ~isempty(ResumeFrom)
%             %Get just the folder path in case user specified a file instead
%             [ResumePath, ResumeName, ResumeExt] = parseFileName(ResumeFrom, 'ignorefilecheck');
%             if isempty(ResumeExt)
%                 ResumePath = fullfile(ResumePath, ResumeName, filesep);
%             end
%             if ~exist(ResumePath, 'dir')
%                 error('%s: Could not find folder to resume from at "%s"', mfilename, ResumeFrom);
%             end
%             showStatus(sprintf('Resuming from %s ...', ResumePath), StatusHandle); 
% 
%             %If *Raw.csv files exist, copy it to temp dir and resume.
%             RawFileStruct = dir(fullfile(ResumePath, '*Raw.csv'));
%             if isempty(RawFileStruct)
%                 error('%s: Could not find *Raw.csv files required for resuming at "%s".', mfilename, ResumePath);
%             end
%             prepTempDir(TempDir);
%             arrayfun(@(x) copyfile(fullfile(ResumePath, x.name), fullfile(TempDir, x.name), 'f'), RawFileStruct);
%             Resume = 'y';
%         end

        %If Resume = 'y', make sure temp dir is not empty
        if strcmpi(Resume, 'y') 
            if ~exist(TempDir, 'dir') || ( exist(TempDir, 'dir') && isempty(dir([TempDir '*Raw.csv'])) )
                warning('%s: Could not resume. Cannot find temp dir with *Raw.csv files at "%s".', mfilename, TempDir); 
                Resume = 'n';
            end
        end

        %Begin initial raw annotation when Resume = 'n' 
        if strcmpi(Resume, 'n')
            %Determine if any temp raw annotation files already exists, which must be deleted
            prepTempDir(TempDir);

            %Determine if you should use segment mode for large files
            if max(SeqRange) == Inf
                SeqCount = getSeqCount(InputFile{f});
            else
                SeqCount = diff(SeqRange);
            end
            if SeqCount > BatchSize
                SegmentMode = 1;
            else
                SegmentMode = 0;
            end

            %Part 1 performs initial annotations in batches
            [~, InFileName, InFileExt] = parseFileName(InputFile{f});
            showStatus(sprintf('Opening %s ...', InFileName), StatusHandle); 
            for j = 1:BatchSize:SeqCount
                %Determine the seq range
                SeqRange = [j j+BatchSize-1];
                if SeqRange(end) > SeqCount
                    SeqRange(end) = SeqCount; 
                end
                if SeqRange(1) < 1
                    SeqRange(1) = 1;
                end
                
                %Correct chain from HL to H if dealing with fasta/q
                if ~isempty(InFileExt) && ismember(lower(InFileExt), {'.fa', '.fasta', '.fastq'})
                    if strcmpi(Chain, 'HL')
                        fprintf('Warning: Only delimited files can do H+L chains. Defaulting to H.\n'); 
                        Chain = 'H';
                    end
                end

                %Open the file and get the sequences within range
                showStatus(sprintf('Processing sequences %d to %d (out of %d) ...', SeqRange(1), SeqRange(end), SeqCount), StatusHandle);
                [VDJdata, VDJheader] = convertInput2VDJdata(InputFile{f}, 'Chain', Chain, 'SeqRange', SeqRange);
                Map = getVDJmapper(VDJheader);

                %Check input sequence for bad characters
                showStatus('Fixing input sequences', StatusHandle);
                [VDJdata, BadIdx] = fixInputSeq(VDJdata, Map);
                if max(BadIdx) ~= 0
%                     error('temp error');
                    saveSeqData(ErrFileName, VDJdata(BadIdx, :), VDJheader, 'append');
                    VDJdata(BadIdx, :) = [];
                end

                %If nothing is left, might be due to wrong delimiter choice
                if isempty(VDJdata)
                    warning('%s: No sequences. Recheck delimiter of input file, which is now set as "%s".', mfilename, Delimiter);
                    continue
                end 

                %Find potential CDR3 start and end locations using V and J gene seed
                %alignment. Do this here, and not when doing VDJ alignment, because users
                %might have complement sequences which must be flipped.
                showStatus('Determining sequence direction and CDR3 areas ...', StatusHandle)
%                 VDJdata = seedCDR3position(VDJdata, Map, DB, 'V', 80,  2, CheckSeqDir);
%                 VDJdata = seedCDR3position(VDJdata, Map, DB, 'J',  3, 14, 'n');
%                 VDJdata = seedCDR3position(VDJdata, Map, DB, 'Vk,Vl', 80,  2, CheckSeqDir);
%                 VDJdata = seedCDR3position(VDJdata, Map, DB, 'Jk,Jl',  3, 14, 'n');

                %Search for initial VDJ alignment matches
                showStatus('Finding initial-guess V(D)J annotations ...', StatusHandle)
                [VDJdata, BadIdx] = findVDJmatch(VDJdata, Map, DB, 'Update', 'Y');
                if max(BadIdx) ~= 0
%                     error('temp error');
                    saveSeqData(ErrFileName, VDJdata(BadIdx, :), VDJheader, 'append');
                    VDJdata(BadIdx, :) = [];
                end

                %Search for initial VJ alignment matches
                [VDJdata, BadIdx] = findVJmatch(VDJdata, Map, DB, 'Update', 'Y');
                if max(BadIdx) ~= 0
%                     error('temp error');
                    saveSeqData(ErrFileName, VDJdata(BadIdx, :), VDJheader, 'append');
                    VDJdata(BadIdx, :) = [];
                end

                %Fix insertion/deletion in V framework
                showStatus('Fixing indels in V genes ...', StatusHandle);
                VDJdata = fixGeneIndel(VDJdata, Map, DB);

                %Remove pseudogenes from degenerate annotations containing functional ones.
                showStatus('Accepting F genes over ORF/P ...', StatusHandle);
                VDJdata = fixDegenVDJ(VDJdata, Map, DB);

                %Insure that V and J segments cover the CDR3 region.
                showStatus('Anchoring 104C and 118W/F ...', StatusHandle);
                VDJdata = constrainGeneVJ(VDJdata, Map, DB);
                
                %Send all Non-functional or incomplete annotation to Err
                showStatus('Moving Nonprod/Invalid/Incomplete Seq to Err file ...', StatusHandle)
                VDJdata = labelSeqQuality(VDJdata, Map, 0.4);
                FunctLoc = [Map.hFunct Map.lFunct];
                FunctLoc(FunctLoc == 0) = [];
                BadLoc = any(cellfun(@(x) contains(x, {'N', 'I', 'M'}), VDJdata(:, FunctLoc)), 2);
                if any(BadLoc)
%                     error('temp error');
                    saveSeqData(ErrFileName, VDJdata(BadLoc, :), VDJheader, 'append');
                    VDJdata(BadLoc, :) = [];
                end
                
                %Finish scheme if annotonly
                if strcmpi(SkipLineage, 'y')
                    VDJdata = findCDR1(VDJdata, Map, DB);
                    VDJdata = findCDR2(VDJdata, Map, DB);
                    VDJdata = buildVDJalignment(VDJdata, Map, DB);
                end

                %Save remaining sequences to temp dir
                if SegmentMode %If sequence file is segmented by CDR3 length
                    CDR3Loc = [Map.hCDR3(1) Map.lCDR3(1)];
                    CDR3Loc(CDR3Loc == 0) = [];
                    CDR3Len = zeros(size(VDJdata, 1), length(CDR3Loc));
                    for k = 1:size(VDJdata, 1)
                        for q = 1:length(CDR3Loc)
                            CDR3Len(k, q) = length(VDJdata{k, CDR3Loc(q)});
                        end
                    end

                    %Save a separate file per unique CDR3H-L combo
                    UnqCDR3Len = unique(CDR3Len, 'rows');
                    for k = 1:size(UnqCDR3Len, 1)
                        Idx = ones(size(CDR3Len, 1), 1, 'logical');
                        for q = 1:length(CDR3Loc)
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

                %If file is small, just save in one file.    
                else
                    SaveName = [TempDir 'Raw.csv'];
                    saveSeqData(SaveName, VDJdata, VDJheader, 'append');
                end
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
            [~, RawFileNameOnly] = parseFileName(RawFileName);
            SrcFile = RawFileName;
            DstFile = fullfile(OutPath, RawFileNameOnly);
            try
                combineSeqData(RawFileList, SrcFile);
                movefile(SrcFile, DstFile);
                delete(RawFileList{:});
            catch
                warning('%s: Could not move Raw files from "%s" to destination "%s".', mfilename, TempDir, DstFile);
            end
        elseif length(RawFileList) == 1
            [~, RawFileNameOnly] = parseFileName(RawFileName);
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