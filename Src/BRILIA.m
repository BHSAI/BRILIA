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
Version = '3.1.3';
varargout = cell(1, nargout);
HasShownCredit = false;

%--------------------------------------------------------------------------
%For running in matlab, make sure BRILIA paths are added correctly
if ~isdeployed
    CurPaths = regexp(path, pathsep, 'split')';
    PathParts = regexp(mfilename('fullpath'), filesep, 'split');
    if ~isempty(PathParts) && isempty(PathParts{1})
        PathParts{1} = filesep;
    end
    MainPath = fullfile(PathParts{1:end-2});
    if ~any(strcmp(CurPaths, MainPath))
        if strcmpi(input('Add BRILIA path to matlab? y or n: ', 's'), 'y')
            fprintf('Adding BRILIA path to Matlab.\n');
            addpath(genpath(MainPath));
        else
            return;
        end
    end
end

%--------------------------------------------------------------------------
%Handle various inputs from matlab or OS command lines
P = inputParser;
addParameter(P, 'InputFile', '', @(x) ischar(x) || iscell(x) || isempty(x));
addParameter(P, 'Chain', 'H', @(x) ismember({upper(x)}, {'H', 'L', 'HL'}));
addParameter(P, 'CheckSeqDir', 'y', @(x) ischar(x) && ismember(lower(x), {'y', 'n'}));
addParameter(P, 'Ddirection', 'all', @(x) ischar(x) && ismember(lower(x), {'all', 'fwd', 'rev', ''}));
addParameter(P, 'DevPerc', 5, @(x) isnumeric(x) && (x>=0) && (x<=100));
addParameter(P, 'Delimiter', '', @(x) ischar(x) && ismember(x, {';', ',', '\t', ''}));
addParameter(P, 'FileType', '', @ischar); %Will make input reader determine file type
addParameter(P, 'NumProc', 'max', @(x) ischar(x) || isnumeric(x));
addParameter(P, 'SettingFile', '', @ischar);
addParameter(P, 'Species', '', @(x) ischar(x) && any(contains(getGeneDatabase('getlist'), x, 'ignorecase', true)));
addParameter(P, 'Strain', 'all', @ischar);
addParameter(P, 'StatusHandle', [], @(x) ishandle(x) || isempty(x) || strcmpi(class(x), 'matlab.ui.control.UIControl'));
addParameter(P, 'Vfunction', 'all', @(x) ischar(x) && min(ismember(regexpi(lower(x), ',', 'split'), {'all', 'f', 'p', 'orf', ''}))==1);
addParameter(P, 'SeqRange', [1,Inf], @(x) isnumeric(x) || ischar(x));
addParameter(P, 'Resume', 'n', @(x) ischar(x) && ismember(lower(x), {'y', 'n'})); %Resumes from known raw file
addParameter(P, 'ResumeFrom', '', @(x) ischar(x)); %File directory storing the *Raw.csv file(s).
addParameter(P, 'OutputFile', [], @(x) ischar(x) || iscell(x) || isempty(x));
addParameter(P, 'BatchSize', 30000, @(x) isnumeric(x) && x >= 1);
addParameter(P, 'AnnotOnly', 'n', @(x) ischar(x) && ismember(lower(x), {'y', 'n'}));
addParameter(P, 'AutoExit', 'y', @(x) ischar(x) && ismember(lower(x), {'y', 'n'}));
addParameter(P, 'SuppressIntro', 'n', @(x) ischar(x) && ismember(lower(x), {'y', 'n'}));

while true
    if ~HasShownCredit && nargin == 0 
        showCredits('bhsai, imgt');
        HasShownCredit = true;
        fprintf('\nType the inputs for BRILIA, ''exit'', or ''help''.\nWarning: File paths with spaces must be wrapped in quotes. Ex: "C:\\Temp Dir\\"\n');
    end
    
    if isempty(varargin)
        Input = input('BRILIA> ', 's');
        if strcmpi(Input, 'exit'); return; end
        QuoteLoc = regexp(Input, '"');
        if mod(length(QuoteLoc), 2) ~= 0
            if nargin == 0
                fprintf('Uneven number of quotes, ".\n');
            else
                error('%s: Uneven number of quotes, ".', mfilename);
            end
        end
        for q = 1:2:length(QuoteLoc)
            Input(QuoteLoc(q):QuoteLoc(q+1)) = strrep(Input(QuoteLoc(q):QuoteLoc(q+1)), ' ', '~');
        end
        CellInput = regexp(Input, '\s+|,', 'split'); 
        varargin = strrep(strrep(CellInput, '~', ' '), '"', '');
    end
    varargin = cleanCommandLineInput(varargin{:});

    %Special parsing of first input
    if ~isempty(varargin) && ischar(varargin{1})
        %Getting version
        if ismember(lower(varargin{1}), {'getversion', 'version'})
            if nargout == 0
                fprintf('BRILIA VERSION %s\n\n', Version);
            else
                varargout{1} = Version;
            end
            if nargin == 0
                varargin = [];
                continue
            else
                return
            end
        end
        
        %Getting help for this
        if ismember(lower(varargin{1}), {'help'})
            showHelp('BRILIA');
            if nargin == 0
                varargin = [];
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
                if nargin == 0
                    varargin = [];
                    continue
                else
                    return
                end
            catch ME
                disp(ME);
            end
        end

        %Check it's a file, correct varargin for parm-value parsing later.
        if exist(varargin{1}, 'file') || exist(fullfile(pwd, varargin{1}), 'file')
            varargin = ['InputFile' varargin];
        end
    end

    try
        [Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInput(P, varargin{:});
    catch ME
        disp(ME);
        varargin = [];
        continue
    end
    if ReturnThis
       varargout = {Ps, Pu, ExpPs, ExpPu};
       return;
    end

    %Override defaults with what is in the SettingFile
    if ~isempty(Ps.SettingFile)
        Ps = readSettingFile(Ps.SettingFile, Ps);
    end

    BatchSize = Ps.BatchSize;
    Chain = Ps.Chain;
    CheckSeqDir = Ps.CheckSeqDir;
    Ddirection = Ps.Ddirection;
    Delimiter = Ps.Delimiter;
    DevPerc = Ps.DevPerc;
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
    SuppressIntro = Ps.SuppressIntro;

    if ~HasShownCredit && ~SuppressIntro
        showCredits('bhsai, imgt');
        HasShownCredit = true;
    end
    
    %--------------------------------------------------------------------------
    %Check the input and output files

    %Get the full file names of input sequences
    if isempty(InputFile) %Ask user to choose
        InputFile = openFileDialog('*.fa*;*.*sv', 'Select the input sequence files', 'multiselect', 'on');
        if isempty(InputFile)
            if nargin == 0
                fprintf('No file was selected.\n\n');
                varargin = [];
                continue
            else
                return
            end
        end
    elseif ischar(InputFile) %Store single file as cell too
        [InputFilePath, InFileName, ~] = parseFileName(InputFile, 'ignorefilecheck'); %Ignore file check for now, as that's done next.
        InputFile = {[InputFilePath InFileName]}; %Ensure file path is always there
    end

    %Make sure same number of output files as input files are specified
    if ~isempty(OutputFile)
        if ischar(OutputFile)
            OutputFile = {OutputFile};
        end
        if length(OutputFile) ~= length(InputFile)
            error('%s: Mismatched number of input (%d) and output files (%d).', mfilename, length(InputFile), length(OutputFile));
        end
    end

    %Check to make sure input files exists
    DelFiles = zeros(1, length(InputFile), 'logical');
    for f = 1:length(InputFile)
        if ~exist(InputFile{f}, 'file')
            DelFiles(f) = 1;
            continue;
        end
        [~, FileName, FileExt] = parseFileName(InputFile{f});
        if ~ismember(lower(FileExt), {'.fa', '.fasta', '.fastq', '.csv', '.tsv', '.ssv', '.txt'})
            warning('%s: Unfamiliar file type (%s) for file (%s). \n  -> Can still process if the "FileType" and "Delimiter" parameters are set.', mfilename, FileExt, FileName);
        end
    end
    InputFile(DelFiles) = [];
    if isempty(InputFile) %No valid files to process
        if nargin ~= 0
            error('%s: No valid input files were provided.', mfilename);
        else
            fprintf('No valid input files were provided.\n');
            varargin = [];
            continue
        end
    end

    %If no output files exists, assign one
    if isempty(OutputFile)
        OutputFile = cell(size(InputFile));
        for f = 1:length(OutputFile)
            [OutputFilePath, OutputFileName, ~] = parseFileName(InputFile{f});
            DotLoc = find(OutputFileName == '.');
            SaveFilePath = [OutputFilePath OutputFileName(1:DotLoc(end)-1) filesep];
            OutputFileName = [OutputFileName(1:DotLoc(end)-1) '.BRILIAv' Version(1) '.csv'];
            OutputFile{f} = [SaveFilePath OutputFileName];
        end
    else
        OutputFile(DelFiles) = [];
    end

    %==========================================================================
    %BRILIA processing begins 

    %Load databases and filter reference genes according to specifications
    if isempty(varargin) %User did not specify species and chain, so ask user.
        [~, ~, FileExt] = parseFileName(InputFile{1});

        %Select the straing
        if ~isempty(FileExt) && ismember(lower(FileExt), {'.fa', '.fasta', '.fastq'})
            ChainList = {'H', 'L'};
            fprintf('Note: only delimited files can do H+L chains.\n'); 
        else
            ChainList = {'H', 'L', 'HL'};
        end
        fprintf('What IG chain is it?\n');
        dispList(ChainList);
        Attempt = 0;
        while true
            Selection = input('Select option: ', 's');
            if isempty(Selection)
                Selection = '1';
            end
            try 
                Selection = round(str2double(Selection));
                if Selection > 0 && Selection <= length(ChainList)
                    Chain = ChainList{Selection};
                    break;
                end
            catch
            end
            Attempt = Attempt + 1;
            if Attempt >= 5
                error('%s: Did not choose correct option.', mfilename);
            end
        end

        %Set these to empty to make the next functions ask the user.
        Species = '';
        Strain = ''; 
        Ddirection = ''; 
        Vfunction = '';
    end
    if strcmpi(SuppressIntro, 'y')
        DB = getGeneDatabase(Species, 'suppress');
    else
        DB = getGeneDatabase(Species);
    end
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
        [OutputFilePath, OutputFileName, ~] = parseFileName(OutputFile{f}, 'ignorefilecheck');
        DotLoc = find(OutputFileName == '.');
        ErrFileName = [OutputFileName(1:DotLoc(end)) 'Err.csv'];
        RawFileName = [OutputFileName(1:DotLoc(end)) 'Raw.csv'];
        TempDir = [OutputFilePath 'Temp' filesep];

        %If using ResumeFrom, try to move Raw.csv into the TempDir.
        if ~isempty(ResumeFrom)
            %Get just the folder path in case user specified a file instead
            [ResumePath, ResumeName, ResumeExt] = parseFileName(ResumeFrom);
            if isempty(ResumePath)
                error('%s: Could not find folder to resume from. \n  "%s"', mfilename, ResumePath);
            end
            if isempty(ResumeExt)
                ResumePath = cat(2, ResumePath, ResumeName, filesep);
            end
            showStatus(sprintf('Resuming from %s ...', ResumePath), StatusHandle); 

            %If *Raw.csv files exist, copy it to temp dir and resume.
            RawFileStruct = dir([ResumePath '*Raw.csv']); 
            if isempty(RawFileStruct)
                error('%s: Could not find Raw file needed for resuming from at "%s"', mfilename, ResumePath);
            end

            %Copy over raw files into the TempDir
            prepTempDir(TempDir);
            for q = 1:length(RawFileStruct)
                copyfile([ResumePath RawFileStruct(q).name], [TempDir RawFileStruct(q).name], 'f');
            end
            Resume = 'y';
        end

        %If Resume = 'y', make sure temp dir is not empty
        if strcmpi(Resume, 'y') 
            if ~exist(TempDir, 'dir') || (exist(TempDir, 'dir') && isempty(dir([TempDir '*Raw.csv'])))
                warning('%s: Could not resume. Cannot find temp dir with *Raw.csv files at "%s".', mfilename, TempDir); 
                Resume = 'n';
            end
        end

        %Begin initial raw annotation of Resume = 'n' 
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
            [~, FileNameOnly, ~] = parseFileName(InputFile{f});
            showStatus(sprintf('Opening %s ...', FileNameOnly), StatusHandle); 
            for j = 1:BatchSize:SeqCount
                %Determine the seq range
                SeqRange = [j j+BatchSize-1];
                if SeqRange(end) > SeqCount
                    SeqRange(end) = SeqCount; 
                end
                if SeqRange(1) < 1
                    SeqRange(1) = 1;
                end

                %Open the file and get the sequences within range
                showStatus(sprintf('Processing sequences %d to %d (total = %d) ...', SeqRange(1), SeqRange(end), SeqCount), StatusHandle);
                [VDJdata, VDJheader] = convertInput2VDJdata(InputFile{f}, 'FileType', FileType, 'Delimiter', Delimiter, 'Chain', Chain, 'SeqRange', SeqRange);
                Map = getVDJmapper(VDJheader);

                %Check input sequence for bad characters
                showStatus('Fixing input sequences', StatusHandle);
                [VDJdata, BadIdx] = fixInputSeq(VDJdata, Map);
                if max(BadIdx) ~= 0
                    saveSeqData([TempDir ErrFileName], VDJdata(BadIdx, :), VDJheader, 'append');
                    VDJdata(BadIdx, :) = [];
                end

                %If nothing is left, might be due to wrong delimiter choice
                if isempty(VDJdata)
                    warning('%s: No sequences. Might have incorrect delimiter. Currently using %s.', mfilename, Delimiter');
                    continue; 
                end 

                %Find potential CDR3 start and end locations using V and J gene seed
                %alignment. Do this here, and not when doing VDJ alignment, because users
                %might have complement sequences which must be flipped.
                showStatus('Determining sequence direction and CDR3 areas ...', StatusHandle)
                VDJdata = seedCDR3position(VDJdata, Map, DB, 'V', 80, 2, CheckSeqDir);
                VDJdata = seedCDR3position(VDJdata, Map, DB, 'J', 3, 14, 'n');
                VDJdata = seedCDR3position(VDJdata, Map, DB, 'Vk, Vl', 80, 2, CheckSeqDir);
                VDJdata = seedCDR3position(VDJdata, Map, DB, 'Jk, Jl', 3, 14, 'n');

                %Search for initial VDJ alignment matches
                showStatus('Finding initial-guess V(D)J annotations ...', StatusHandle)
                [VDJdata, BadIdx] = findVDJmatch(VDJdata, Map, DB, 'Update', 'Y');
                if max(BadIdx) ~= 0
                    saveSeqData([TempDir ErrFileName], VDJdata(BadIdx, :), VDJheader, 'append');
                    VDJdata(BadIdx, :) = [];
                end

                %Search for initial VJ alignment matches
                [VDJdata, BadIdx] = findVJmatch(VDJdata, Map, DB, 'Update', 'Y');
                if max(BadIdx) ~= 0
                    saveSeqData([TempDir ErrFileName], VDJdata(BadIdx, :), VDJheader, 'append');
                    VDJdata(BadIdx, :) = [];
                end

                %Send all Non-functional to error
                FunctLoc = [Map.hFunct Map.lFunct];
                FunctLoc(FunctLoc == 0) = [];
                BadIdx = zeros(size(VDJdata, 1), 1, 'logical');
                for w = 1:size(VDJdata, 1)
                    for k = 1:length(FunctLoc)
                        if isempty(VDJdata{w, FunctLoc(k)})
                            BadIdx(w) = 1;
                            break;
                        end
                        if ~isempty(VDJdata{w, FunctLoc(k)}) && VDJdata{w, FunctLoc(k)} ~= 'Y'
                            BadIdx(w) = 1;  
                            break;
                        end
                    end
                end
                if max(BadIdx) ~= 0
                    saveSeqData([TempDir ErrFileName], VDJdata(BadIdx, :), VDJheader, 'append');
                    VDJdata(BadIdx, :) = [];
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
            clear VDJdata;
        end

        %Part 2 does everything AFTER intial annotations
        GrpNumStart = 0 ;
        FileList = dir([TempDir '*Raw.csv']);

        if strcmpi(AnnotOnly, 'n')
            for t = 1:length(FileList)
                %Check if a Final.csv file already exists to skip
                TempRawFileName = FileList(t).name;
                TempOutputFileName = strrep(TempRawFileName, 'Raw.csv', 'Final.csv');
                if exist(TempOutputFileName, 'file')
                    continue
                end

                %Reload the sequence file with same-length CDR3s
                showStatus(sprintf('Processing %s ...', TempRawFileName), StatusHandle);
                [VDJdata, VDJheader] = openSeqData([TempDir TempRawFileName]);
                Map = getVDJmapper(VDJheader);

                %Fix insertion/deletion in V framework
                showStatus('Fixing indels in V genes ...', StatusHandle)
                VDJdata = fixGeneIndel(VDJdata, Map, DB);

                %Remove pseudogenes from degenerate annotations containing functional ones.
                showStatus('Accepting F genes over ORF/P ...', StatusHandle)
                VDJdata = fixDegenVDJ(VDJdata, Map, DB);

                %Insure that V and J segments cover the CDR3 region.
                showStatus('Anchoring 104C and 118W/F ...', StatusHandle)
                VDJdata = constrainGeneVJ(VDJdata, Map, DB);

                %Cluster the data based variable region and hamming dist of DevPerc%.
                showStatus('Clustering by lineage ...', StatusHandle)
                [VDJdata, BadVDJdataT] = clusterGene(VDJdata, Map, DevPerc);
                if size(BadVDJdataT, 1) ~= 0
                    saveSeqData([TempDir ErrFileName], BadVDJdataT, VDJheader, 'append');
                    clear BadVDJdataT;
                end
                if size(VDJdata, 1) == 0; continue; end

                %Renumbering groups
                GrpNums = cell2mat(VDJdata(:, Map.GrpNum)) + GrpNumStart;
                VDJdata(:, Map.GrpNum) = num2cell(GrpNums);
                GrpNumStart = max(GrpNums); 

                %Set all groups to have same annotation and VMDNJ lengths.
                showStatus('Correcting annotations by lineage ...', StatusHandle)
                VDJdata = conformGeneGroup(VDJdata, Map, DB);

                %Get better D match based on location of consensus V J mismatches.
                showStatus('Refining D annotations ...', StatusHandle)
                VDJdata = findBetterD(VDJdata, Map, DB);

                %Trim V, D, J edges and extract better N regions
                showStatus('Trimming N regions ...', StatusHandle)
                VDJdata = trimGeneEdge(VDJdata, Map, DB);

                %Fix obviously incorrect trees.
                showStatus('Rerooting lineage trees ...', StatusHandle)
                VDJdata = removeDupSeq(VDJdata, Map);
                VDJdata = fixTree(VDJdata, Map);

                %Finalize VDJdata details and CDR 1, 2, 3 info
                VDJdata = padtrimSeqGroup(VDJdata, Map, 'grpnum', 'trim', 'Seq'); %will only remove "x" before and after Seq if they all have it. 
                VDJdata = findCDR1(VDJdata, Map, DB);
                VDJdata = findCDR2(VDJdata, Map, DB);
                VDJdata = findCDR3(VDJdata, Map, DB, 'IMGT'); %removes the 104C and 118W from CDR3, and adjusts the CDR3 length to true IMGT length

                %Move non-function sequences to Err file too
                FunctLoc = [Map.hFunct; Map.lFunct];
                FunctLoc(FunctLoc == 0) = [];
                BadIdx = zeros(size(VDJdata, 1), 1, 'logical');
                for w = 1:size(VDJdata, 1)
                    for q = 1:length(FunctLoc)
                        if isempty(VDJdata{w, FunctLoc(q)}) || VDJdata{w, FunctLoc(q)} == 'N'
                            BadIdx(w) = 1;
                            break
                        end
                    end
                end
                if max(BadIdx) ~= 0
                    saveSeqData([TempDir ErrFileName], VDJdata(BadIdx, :), VDJheader, 'append');
                end
                VDJdata(BadIdx, :) = [];

                %Save the functional annotations
                VDJdata = buildVDJalignment(VDJdata, Map, DB); %Adds the alignment information
                saveSeqData([TempDir TempOutputFileName], VDJdata, VDJheader, 'append');
                clear VDJdata BadVDJdata
            end
            %======================================================================
            %Move folder to the correct destination folder

            %Combine final annotations to a single file, then move to destination
            FinalFileStruct = dir([TempDir '*Final.csv']);
            FinalFileList = cell(length(FinalFileStruct), 1);
            for q = 1:length(FinalFileList)
                FinalFileList{q} = [TempDir FinalFileStruct(q).name];
            end
            if exist([TempDir OutputFileName], 'file') %If a real output file was interrupted and left in Temp folder, delete.
                try
                    delete([TempDir OutputFileName]);
                catch
                    warning('%s: Could not delete incomplete final file [ %s ].', mfilename, [TempDir OutputFileName]);
                end
            end
            if length(FinalFileList) > 1
                combineSeqData(FinalFileList, [TempDir OutputFileName]);
                try
                    movefile([TempDir OutputFileName], [OutputFilePath, OutputFileName], 'f');
                    delete(FinalFileList{:});
                catch
                    warning('%s: Could not move Final files from [ %s ] to destination [ %s ].', mfilename, [TempDir OutputFileName], [OutputFilePath OutputFileName]);
                end
            elseif length(FinalFileList) == 1
                try
                    movefile(FinalFileList{1}, [OutputFilePath, OutputFileName], 'f');
                catch
                    warning('%s: Could not move Final files from [ %s ] to destination [ %s ].', mfilename, FinalFileList{1}, [OutputFilePath OutputFileName]);
                end
            end
        end

        %Combine raw annotations to a single file, then move to destination
        RawFileStruct = dir([TempDir '*Raw.csv']); 
        RawFileList = cell(length(RawFileStruct), 1);
        for q = 1:length(RawFileList)
            RawFileList{q} = [TempDir RawFileStruct(q).name];
        end
        if length(RawFileList) > 1
            try
                combineSeqData(RawFileList, [TempDir RawFileName]);
                movefile([TempDir RawFileName], [OutputFilePath, RawFileName]);
                delete(RawFileList{:});
            catch
                warning('%s: Could not move Raw files from [ %s ] to destination [ %s ].', mfilename, TempDir, [OutputFilePath RawFileName]);
            end
        elseif length(RawFileList) == 1
            try
                movefile(RawFileList{1}, [OutputFilePath, RawFileName], 'f');
            catch
                warning('%s: Could not move Raw files from [ %s ] to destination [ %s ].', mfilename, RawFileList{1}, [OutputFilePath RawFileName]);
            end
        end

        %Move the err file out
        if exist([TempDir ErrFileName], 'file')
            try
                movefile([TempDir ErrFileName], [OutputFilePath, ErrFileName], 'f');
            catch
                warning('%s: Could not move Err files from [ %s ] to destination [ %s ].', mfilename, [TempDir ErrFileName], [OutputFilePath ErrFileName]);
            end
        end

        %If isempty TempDir, delete
        prepTempDir(TempDir, 'delete');

        showStatus(sprintf('Finished in %0.1f sec.', toc), StatusHandle);
    end
    varargout{1} = OutputFile;
    
    if strcmpi(AutoExit, 'y'); return; end
    varargin = []; %Empty this to make sure loop asks for next input.
end