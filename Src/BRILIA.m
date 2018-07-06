%                              BRILIA
%(B-cell Repertoire Inductive Lineage and Immunosequence Annotator)
%Written by Donald Lee (dlee@bhsai.org)
%
%After starting BRILIA in the command line envinronment, directly type the
%inputs as shown.
%
%  Opens a GUI version of BRILIA
%    > GUI_BRILIA
%  
%  Opens a GUI version of plotTree function
%    > GUI_plotTree
%
%  Runs annotation
%    > [SequenceFile.fasta/fastq/csv] Param1 Value1 Param2 Value2 ...
%
%  Param-Value are case-insensitive and shown below, ordered by category (C) and relevance. 
%    Categories (C) are: d=database, p=performance, i=input, o=output, l=lineage.
%    * = default value, [] = variable string, and # = variable number.
%
%    C  Param         Value                      Details
%    -  ------------  -------------------------  -----------------------------
%    d  Chain         h* | l | hl                IgG heavy and/or light chain.
%    d  Species       human | mouse | macaque    Database by species.
%                     zebrafish
%    d  Strain        all* | c57bl | balb        Database by strain, only for mouse.
%    d  CheckSeqDir   y* | n                     Check for rev-comp seq?
%    d  Ddirection    all* | fwd | rev           D gene direction.
%    d  Vfunction     all* | f | p | orf         V gene functions.
%    l  DevPerc       # (btwn. 0 and 100)        Lineage cutoff hamming dist, SHM-adjusted (deviation percentage).
%    p  BatchSize     # | 30000*                 Number of sequences to process in batches to prevent memory overload.
%    p  NumProc       max | #                    Set max or # of cores to use.
%    i  InputFile     ''* | [InputFile.*]        Input sequence file (.fasta or .fastaq or .csv). Default is empty and will ask  user to choose a file.
%    i  Delimiter     , | ; | \t                 (will autodetect)
%    i  FileType      fasta | fastq | delimited  (will autodetect)
%    i  Resume        y | n*                     Resume interrupted job? Must have files in the Temp folder and use the same database.
%    i  ResumeFrom    [FolderName]               Resume from Raw.csv files from a prior job's folder containg the Raw.csv. Must use same database.
%    i  SeqRange      # | [#,#] (incl. brackets) The #th sequence, or #th to #th sequences to process. Default, does all.
%    i  SettingFile   [SettingsFile.txt]         Preset setting txt file.
%    o  AnnotOnly     y | n*                     Do only annotation and skip lineage-based correction.                   
%    o  OutputFile    [OutputFile.csv]           Specify if default output file names are unwanted.
%
%  OUTPUT FILES
%    [filename]BRILIAvN.csv     : final annotation file
%    [filename]BRILIAvN.Raw.csv : annotation file w/o lineage correction
%    [filename]BRILIAvN.Err.csv : error or non-productive annotations
%
%  OUTPUT DIR
%    [InputFilePath/InputFileName/]: Output files will be stored in the
%      same location as the input file, but under a folder with the same
%      name as the input file.
%    [InputFilePath/Temp] : temporary folder storing annotation files. This
%      serves as a save point too when BRILIA is interrupted due to
%      external sources (system shutdown). In such case, use the same input
%      command BUT add 'Resume y' option to continue annotation.
%
%  EXAMPLE
%    > SeqFile.fastq SettingFile SettingFile.txt
%    > SeqFile.fasta Species Mouse Strain C57BL Chain H CheckSeqDir Y NumProc 4 SeqRange [1,100]

function varargout = BRILIA(varargin)
Version = '3.2.0';
varargout = cell(1, nargout);
HasShownCredit = false;

%--------------------------------------------------------------------------
%For running in matlab, make sure BRILIA paths are added correctly
if ~isdeployed
    CurPaths = regexp(path, pathsep, 'split')';
    MainPath = fileparts(fileparts(mfilename('fullpath')));
    if ~any(strcmpi(CurPaths, MainPath))
        if strcmpi(input('Add BRILIA path to matlab? y or n: ', 's'), 'y')
            fprintf('Adding BRILIA path to Matlab.\n');
            addpath(genpath(MainPath));
        else
            return
        end
    end
end

%--------------------------------------------------------------------------
%Handle various inputs from matlab or OS command lines
P = inputParser;
addParameter(P, 'InputFile',     '',      @(x) ischar(x) || iscell(x) || isempty(x));
addParameter(P, 'Chain',         'H',     @(x) ismember({upper(x)}, {'H', 'L', 'HL'}));
addParameter(P, 'CheckSeqDir',   'y',     @(x) ischar(x) && ismember(lower(x), {'y', 'n'}));
addParameter(P, 'Ddirection',    'all',   @(x) ischar(x) && ismember(lower(x), {'all', 'fwd', 'rev', ''}));
addParameter(P, 'Delimiter',     '',      @(x) ischar(x) && ismember(x, {';', ',', '\t', ''}));
addParameter(P, 'FileType',      '',      @ischar); %Will make input reader determine file type
addParameter(P, 'NumProc',       'max',   @(x) ischar(x) || isnumeric(x));
addParameter(P, 'SettingFile',   '',      @ischar);
addParameter(P, 'Species',       '',      @(x) ischar(x) && any(contains(getGeneDatabase('getlist'), x, 'ignorecase', true)));
addParameter(P, 'Strain',        'all',   @ischar);
addParameter(P, 'StatusHandle',  [],      @(x) ishandle(x) || isempty(x) || strcmpi(class(x), 'matlab.ui.control.UIControl'));
addParameter(P, 'Vfunction',     'all',   @(x) ischar(x) && min(ismember(regexpi(lower(x), ',', 'split'), {'all', 'f', 'p', 'orf', ''}))==1);
addParameter(P, 'SeqRange',      [1,Inf], @(x) isnumeric(x) || ischar(x));
addParameter(P, 'Resume',        'n',     @(x) ischar(x) && ismember(lower(x), {'y', 'n'})); %Resumes from known raw file
addParameter(P, 'ResumeFrom',    '',      @(x) ischar(x)); %File directory storing the *Raw.csv file(s).
addParameter(P, 'OutputFile',    [],      @(x) ischar(x) || iscell(x) || isempty(x));
addParameter(P, 'BatchSize',     30000,   @(x) isnumeric(x) && x >= 1);
addParameter(P, 'AnnotOnly',     'n',     @(x) ischar(x) && ismember(lower(x), {'y', 'n'}));
addParameter(P, 'AutoExit',      'y',     @(x) ischar(x) && ismember(lower(x), {'y', 'n'}));

%Determine the BRILIA run mode
if nargin == 0
    RunMode = 3;     %Running from EXE or Matlab without inputs
else
    if isdeployed
        RunMode = 2; %Running from EXE with inputs
    else
        RunMode = 1; %Running from matlab with inputs
    end
end

while true
    if RunMode == 3
        if ~HasShownCredit
            showCredits('bhsai', 'imgt');
            HasShownCredit = true;
            fprintf('\nType input commands, ''exit'', or ''help''.  (Wrap files with spaces in quotes, like "C:\\Temp Dir\\)":\n'); %Only show this once
        end
        
        Input = input('BRILIA> ', 's'); %If empty, will ask user to choose file
        if strcmpi(Input, 'exit')
            return
        end
    else
        for k = 1:length(varargin)
            if isnumeric(varargin{k})
                varargin{k} = mat2str(varargin{k});
            else
                varargin{k} = strrep(varargin{k}, ' ', '^'); %Sub spaces with ^ to prevent error parsing
            end
        end
        if nargin > 1
            Input = [sprintf('%s ', varargin{1:end-1}) varargin{end}];
        else
            Input = varargin{1};
        end
    end
        
    %Check to make sure there are paired number of quotes ".
    QuoteLoc = regexp(Input, '"');
    if mod(length(QuoteLoc), 2) > 0
        if RunMode == 3
            fprintf('Error: Uneven number of quotes ( " ) in inputs.\n');
            continue
        else
            error('%s: Uneven number of quotes ( " ) in inputs.\n', mfilename);
        end
    end
        
    %Format string inputs with quotes to proper BRILIA inputs
    for q = 1:2:length(QuoteLoc)
        Input(QuoteLoc(q):QuoteLoc(q+1)) = strrep(Input(QuoteLoc(q):QuoteLoc(q+1)), ' ', '^');
    end
    CellInput = regexp(Input, '\s+|,', 'split'); 
    varargin = strrep(strrep(CellInput, '^', ' '), '"', '');
    varargin = cleanCommandLineInput(varargin{:});

    %Special parsing of first input
    if ~isempty(varargin) && ischar(varargin{1})
        %Getting version
        if ismember(lower(varargin{1}), {'getversion', 'version'})
            if RunMode == 3
                fprintf('BRILIA VERSION %s\n\n', Version);
                continue
            else
                varargout{1} = Version;
                return
            end
        end
        
        %Getting help for this
        if ismember(lower(varargin{1}), {'h', 'i', 'info', 'help', 'showhelp', 'showinfo'})
            if length(varargin) == 1
                showHelp('BRILIA');
            else
                showHelp(varargin{2:end});
            end
            if RunMode == 3
                continue
            else
                return
            end
        end

        %Redirect to BRILIA subfunction
        SubFuncNames = {'pwd', 'cd', 'dir', 'ls', 'plotTree', 'runAnalysis', 'setCores', 'showHelp', 'GUI_BRILIA', 'GUI_plotTree'}; %For security, only accept allowed function calls!
        if ismember(varargin{1}, SubFuncNames)
            try
                FH = str2func(varargin{1});
                FH(varargin{2:end});
            catch ME
                disp(ME);
            end
            if RunMode == 3
                continue
            else
                return
            end
        end

        %Check it's a file, correct varargin for parm-value parsing later.
        if exist(varargin{1}, 'file') || isempty(varargin{1})
            varargin = ['InputFile' varargin];
        end
    end

    %Parse the inputs and return if this was a code-2-code summon
    try
        [Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInput(P, varargin{:});
        if ReturnThis
           varargout = {Ps, Pu, ExpPs, ExpPu};
           return
        end
    catch ME
        fprintf('%s: Error parsing input at line = %d.\n', mfilename, ME.stack(1).line);
        if RunMode == 3
            continue
        else
            return
        end
    end

    %Override defaults with what is in the SettingFile
    if ~isempty(Ps.SettingFile)
        Ps = readSettingFile(Ps.SettingFile, Ps);
    end

    BatchSize = Ps.BatchSize;
    MainChain = Ps.Chain;
    CheckSeqDir = Ps.CheckSeqDir;
    Ddirection = Ps.Ddirection;
    Delimiter = Ps.Delimiter;
    FileType = Ps.FileType;
    InputFile = Ps.InputFile;
    NumProc = Ps.NumProc;
    OutputFile = Ps.OutputFile;
    Resume = Ps.Resume;
    ResumeFrom = Ps.ResumeFrom;
    SeqRange = Ps.SeqRange;
    Species = Ps.Species;
    StatusHandle = Ps.StatusHandle;
    Strain = Ps.Strain;
    Vfunction = Ps.Vfunction;
    AnnotOnly = Ps.AnnotOnly;
    AutoExit = Ps.AutoExit;

    %Show credits AFTER input parsing for runmode 1 or 2
    if ~HasShownCredit
        showCredits('bhsai', 'imgt');
        HasShownCredit = true;
    end
    
    %--------------------------------------------------------------------------
    %Check the input and output files

    %Get the full file names of input sequences
    if isempty(InputFile) %Ask user to choose
        InputFile = openFileDialog('*.fa*;*.*sv', 'Select the input sequence files', 'multiselect', 'on');
        if isempty(InputFile)
            if RunMode == 3
                fprintf('No file was selected.\n\n');
                continue
            else
                return
            end
        end
    elseif ischar(InputFile) %Store single file as cell too
        [InputFilePath, InFileName, ~] = parseFileName(InputFile, 'ignorefilecheck'); %Ignore file check for now, as that's done next.
        InputFile = {fullfile(InputFilePath, InFileName)}; %Ensure file path is always there
    end

    %Make sure # of output files = # of input files
    if ~isempty(OutputFile)
        if ischar(OutputFile)
            OutputFile = {OutputFile};
        end
        if length(OutputFile) ~= length(InputFile)
            error('%s: Mismatched number of input (%d) and output files (%d).', mfilename, length(InputFile), length(OutputFile));
        end
    end

    %Delete input files that do not exist
    DelFiles = zeros(1, length(InputFile), 'logical');
    for f = 1:length(InputFile)
        if ~exist(InputFile{f}, 'file')
            DelFiles(f) = 1;
            continue
        end
        [~, FileName, FileExt] = parseFileName(InputFile{f});
        if ~ismember(lower(FileExt), {'.fa', '.fasta', '.fastq', '.csv', '.tsv', '.ssv', '.txt'})
            warning('%s: Unfamiliar file type "%s" for file "%s". \n  -> Can still process if the "FileType" and "Delimiter" parameters are set.', mfilename, FileExt, FileName);
        end
    end
    if any(DelFiles)
        fprintf('%s: Could not find the following input files:\n', mfilename);
        fprintf('  %s\n', InputFile{DelFiles});
        InputFile(DelFiles) = [];
    end
    
    if isempty(InputFile)
        if RunMode == 3
            fprintf('No valid input files were provided.\n');
            continue
        else
            error('%s: No valid input files were provided.', mfilename);
        end
    end

    %If no output files exists, assign one
    if isempty(OutputFile)
        OutputFile = cell(size(InputFile));
        for f = 1:length(OutputFile)
            [OutPath, ~, ~, OutFilePre] = parseFileName(InputFile{f});
            OutputFile{f} = fullfile(OutPath, OutFilePre, [OutFilePre '.BRILIAv' Version(1) '.csv']);
        end
    else
        OutputFile(DelFiles) = [];
    end

    %==========================================================================
    %BRILIA processing begins 

    %See if the user specified everything when using RunMode = 3
    if RunMode == 3 
        %Selecting the H or L chain
        if ~strcmpi('Chain', varargin) 
            fprintf('What IG chain is it?\n');
            ChainList = {'H', 'L', 'HL'};
            dispList(ChainList);
            Attempt = 0;
            while true
                try 
                    Selection = round(input('Select option (no selection = 1): '));
                    if isempty(Selection)
                        Selection = 1;
                        break
                    elseif Selection >= 1 && Selection <= length(ChainList)
                        break
                    end
                    Attempt = Attempt + 1;
                catch
                    Attempt = Attempt + 1;
                end
                if Attempt >= 5
                    fprintf('Error: Did not choose correct option.\n');
                    break
                end
            end
            if Attempt >= 5
                continue %return to beginning of input because of failure
            end
            MainChain = ChainList{Selection};
        end
        
        if ~strcmpi('Species', varargin) 
            Species = '';
        end
        if ~strcmpi('Strain', varargin)
            Strain = '';
        end
        if ~strcmpi('Ddirection', varargin) 
            Ddirection = '';
        end
        if ~strcmpi('Vfunction', varargin) 
            Vfunction = '';
        end
    end
    
    %Load databases and filter reference genes according to specifications
    DB = getGeneDatabase(Species);
    DB = filterGeneDatabase(DB, 'Strain', Strain, 'Ddirection', Ddirection, 'Vfunction', Vfunction);

    %Set the number of processors to use
    showStatus('Setting up parallel computing ...', StatusHandle);
    [~, NumWorkers] = setCores(NumProc);
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

        %If using ResumeFrom, try to move Raw.csv into the TempDir.
        if ~isempty(ResumeFrom)
            %Get just the folder path in case user specified a file instead
            [ResumePath, ResumeName, ResumeExt] = parseFileName(ResumeFrom, 'ignorefilecheck');
            if isempty(ResumeExt)
                ResumePath = fullfile(ResumePath, ResumeName, filesep);
            end
            if ~exist(ResumePath, 'dir')
                error('%s: Could not find folder to resume from at "%s"', mfilename, ResumeFrom);
            end
            showStatus(sprintf('Resuming from %s ...', ResumePath), StatusHandle); 

            %If *Raw.csv files exist, copy it to temp dir and resume.
            RawFileStruct = dir(fullfile(ResumePath, '*Raw.csv'));
            if isempty(RawFileStruct)
                error('%s: Could not find *Raw.csv files required for resuming at "%s".', mfilename, ResumePath);
            end
            prepTempDir(TempDir);
            arrayfun(@(x) copyfile(fullfile(ResumePath, x.name), fullfile(TempDir, x.name), 'f'), RawFileStruct);
            Resume = 'y';
        end

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
                Chain = MainChain;
                if ~isempty(InFileExt) && ismember(lower(FileExt), {'.fa', '.fasta', '.fastq'})
                    if strcmpi(Chain, 'HL')
                        fprintf('Warning: Only delimited files can do H+L chains. Defaulting to H.\n'); 
                        Chain = 'H';
                    end
                end

                %Open the file and get the sequences within range
                showStatus(sprintf('Processing sequences %d to %d (out of %d) ...', SeqRange(1), SeqRange(end), SeqCount), StatusHandle);
                [VDJdata, VDJheader] = convertInput2VDJdata(InputFile{f}, 'FileType', FileType, 'Delimiter', Delimiter, 'Chain', Chain, 'SeqRange', SeqRange);
                Map = getVDJmapper(VDJheader);

                %Check input sequence for bad characters
                showStatus('Fixing input sequences', StatusHandle);
                [VDJdata, BadIdx] = fixInputSeq(VDJdata, Map);
                if max(BadIdx) ~= 0
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
                VDJdata = seedCDR3position(VDJdata, Map, DB, 'V', 80,  2, CheckSeqDir);
                VDJdata = seedCDR3position(VDJdata, Map, DB, 'J',  3, 14, 'n');
                VDJdata = seedCDR3position(VDJdata, Map, DB, 'Vk,Vl', 80,  2, CheckSeqDir);
                VDJdata = seedCDR3position(VDJdata, Map, DB, 'Jk,Jl',  3, 14, 'n');

                %Search for initial VDJ alignment matches
                showStatus('Finding initial-guess V(D)J annotations ...', StatusHandle)
                [VDJdata, BadIdx] = findVDJmatch(VDJdata, Map, DB, 'Update', 'Y');
                if max(BadIdx) ~= 0
                    saveSeqData(ErrFileName, VDJdata(BadIdx, :), VDJheader, 'append');
                    VDJdata(BadIdx, :) = [];
                end

                %Search for initial VJ alignment matches
                [VDJdata, BadIdx] = findVJmatch(VDJdata, Map, DB, 'Update', 'Y');
                if max(BadIdx) ~= 0
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
                    saveSeqData(ErrFileName, VDJdata(BadLoc, :), VDJheader, 'append');
                    VDJdata(BadLoc, :) = [];
                end
                
                %Finish scheme if annotonly
                if strcmpi(AnnotOnly, 'y')
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
        GrpNumStart = 1 ;
        FileList = dir([TempDir '*Raw.csv']);

        if strcmpi(AnnotOnly, 'n')
            for t = 1:length(FileList)
                %Check if a Final.csv file already exists to skip
                TempRawFileName = FileList(t).name;
                TempOutFileName = strrep(TempRawFileName, 'Raw.csv', 'Final.csv');
                if exist(fullfile(TempDir, TempOutFileName), 'file'); continue; end
                
                %Reload the sequence file with same-length CDR3s
                showStatus(sprintf('Processing %s ...', TempRawFileName), StatusHandle);
                [VDJdata, VDJheader] = openSeqData(fullfile(TempDir, TempRawFileName));
                Map = getVDJmapper(VDJheader);
                if size(VDJdata, 1) < 1; continue; end
                
                %Cluster the data based variable region and hamming dist of DevPerc%.
                showStatus('Clustering by lineage ...', StatusHandle)
                VDJdata = clusterGene(VDJdata, Map);

                %Renumbering groups since they're done in batches
                GrpNums = cell2mat(VDJdata(:, Map.GrpNum)) + GrpNumStart - 1;
                VDJdata(:, Map.GrpNum) = num2cell(GrpNums);
                GrpNumStart = max(GrpNums) + 1; 

                %Set all groups to have same annotation and VMDNJ lengths.
                showStatus('Correcting annotations by lineage ...', StatusHandle)
                VDJdata = conformGeneGroup(VDJdata, Map, DB);

                %Get better D match based on location of consensus V J mismatches.
                showStatus('Refining D annotations ...', StatusHandle)
                VDJdata = findBetterD(VDJdata, Map, DB);

                %Trim V, D, J edges and extract better N regions
                showStatus('Trimming N regions ...', StatusHandle)
                VDJdata = trimGeneEdge(VDJdata, Map, DB);

                %Finalize VDJdata details and CDR 1, 2, 3 info
                VDJdata = padtrimSeqGroup(VDJdata, Map, 'grpnum', 'trim', 'Seq'); %will only remove "x" before and after Seq if they all have it. 
                VDJdata = findCDR1(VDJdata, Map, DB);
                VDJdata = findCDR2(VDJdata, Map, DB);
                VDJdata = findCDR3(VDJdata, Map, DB, 'IMGT'); %removes the 104C and 118W from CDR3, and adjusts the CDR3 length to true IMGT length
                VDJdata = buildVDJalignment(VDJdata, Map, DB);
                
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
    
    if RunMode == 3 && strcmpi(AutoExit, 'n')
        continue 
    else 
        return
    end
end