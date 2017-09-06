%BRILIA is the main code that will run the annotation of immunosequencing
%data. BRILIA stands for B-cell Repertoire Inductive Lineage and
%Immunosequence Annotator, where "immunosequence" is simply referring to a
%DNA/RNA sequences of a B cell receptor. BRILIA works better if the
%sequences are obtained from a single repertoire of B cells from the same
%host, and if low-quality regions are trimmed away or replaced with "N" or
%"X" wildcard nucleotides.
%
%  BRILIA('getversion')
%
%  BRILIA(InputFile, Param, Value, ...)
%
%  BRILIA(FunctionName, ...)
%
%  INPUT
%    InputFile: char or cell name(s) of the sequence file(s) to process
%    FunctionName: a BRILIA internal function that can be summoned, in
%    which following inputs are specific to this internal function. This
%    allows users to summon plotTree function via BRILIA's binary file. See
%    NOTE.
%
%    Param, Value pairs are for specific setting as follows
%      Param             Value                  This param-value sets:
%      ---------------   ---------------------  -----------------------------
%      'SettingFile'     [SettingsFile.txt]     All settings specified by
%                                                 a txt file
%      'Species'         'human' 'mouse' etc    VDJ database by species
%      'Strain'          'all' 'C57BL' etc      VDJ database by strain
%      'Ddirection'      'all' 'fwd' 'rev'      Allowed D gene direction
%      'Vfunction'       'all' 'f' 'p' 'orf'    Allowed V gene functions
%      'DevPerc'         0 <= N <= 100          Clustering cutoff dist by
%                                                 of seq length
%      'FileType'        'fasta', 'fastq',      Speficy if input file type 
%                        'delimited'              unobvious (eg, seq.txt)
%      'Delimiter'       ';' ', ' '\t' ''       Delimiter type. Empty means
%                                                 autodetect. Unused if
%                                                 loading fasta or fastq.
%      'CheckSeqDir'     'y' 'n'                Whether or not to check 
%                                                 for complement seq
%      'NumProc'         'max' or N             Number of processors to
%                                                 use. Default 'max'.
%      'StatusHandle'    [text handle obj]      GUI handle used to update
%                                                 the status text.
%      'SeqRange'        N or [N(1) N(2)]       Nth or N(1) to N(2)
%                                                 sequence to process.
%      'Resume'          'y' or 'n'             'y' = start from beginning
%                                               'n' = resumes from an
%                                                 interrupted job based on
%                                                 what is in the Temp dir.
%      'ResumeFrom'      [FolderName]           Directory storing the
%                                                 Raw.csv file to resume
%                                                 from for the lineage
%                                                 correction. Skips
%                                                 annotation, but must use
%                                                 same database.
%      'OutputFile'      [OutputFile.csv]       Outputfile name. See note
%                                                 below about outputfile
%                                                 handling.
%
%  OUTPUT
%    Version: Version string in X.Y.Z
%    P: default inputs of BRILIA returned as a structure
%    SaveFileName: the final annotation file names
%    [filename].Raw.csv: annotation only file, no lineage correction. Used for
%    resuming from ('ResumeFrom') annotation file to lineage-corrected
%    annotation file.
%    [filename].Err.csv: sequences that could not be annotated fully or are
%    unproductive.
%    [filename].csv: final annotation file
%
%  NOTE
%    All OutputFile is a comma-delimited ".csv" file.
%
%    A Temp folder will ALWAYS be created in the same folder where the
%    InputFile is. This folder stores the initial VDJ annotation files,
%    which is used as a save point for conducting lineage correction in
%    case there was an interruption. Without this folder, BRILIA cannot
%    "resume" from annotation and must start from the beginning, which
%    takes a while, especially if one just want to change the lineage
%    parameter.
%
%    If NO OutputFile is given, then it will save the file into the input
%    file directory, subfolder "BRILIA", under the following filename,
%    inputfilename.BRILIAvX.[XXX].csv. This file is
%    semicolon delimited.
%    
%    If OutputFile WITH the full path is given, then it will create the
%    folder directory of the output file and save everything to this file
%    when done. Note that an [OutputFile].Raw.csv file will be created,
%    but the final annotations is [OutputFile].csv and NOT
%    [OutputFile].Final.csv.
%    
%    If OutputFile WITHOUT the full path is given, then it will create a
%    subfolder called "BRILIA" where the input file is, and save the
%    outputfile there, along with the [OutputFile].Raw.csv file and
%    [OutputFile].csv.
%
%    Users can summon another function within BRILIA by inputting the
%    function name, which follows a lower-case camel naming convention such
%    as plotTree. Example: BRILIA('plotTree', BRILIA_OuputFile, 'DPI', 600)
%    will summon the plotTree function, open BRILIA_OutputFile, and save
%    images as 600 DPI to a default output folder specified by the plotTree
%    function. This capability was added to reduce the file size of BRILIA
%    as codes are redundant.
%
%  EXAMPLE (MATLAB)
%    >> BRILIA([], 'SettingFile', 'SettingFile.txt');
%
%    >> BRILIA('SeqFile.fasta', 'Species', 'Mouse', 'Strain', 'C57BL',
%      'Chain', 'H', 'CheckSeqDir', 'n', 'NumProc', '4');
%
%  EXAMPLE (command line)
%    > BRILIA SeqFile.fasta SettingFile SettingFile.txt
%
%    > BRILIA SeqFile.fasta Species Mouse Strain C57BL Chain H CheckSeqDir
%      N NumProc 4
%
%  Written by Donald Lee, dlee@bhsai.org

function varargout = BRILIA(varargin)
Version = '3.0.9';

%--------------------------------------------------------------------------
%For running in matlab, make sure BRILIA paths are added correctly
if ~isdeployed
    CurPaths = regexp(path, pathsep, 'split')';
    MainPath = mfilename('fullpath');
    SlashLoc = regexp(MainPath, filesep);
    MainPath = MainPath(1:SlashLoc(end)-1);
    HavePath = 0;
    for p = 1:length(CurPaths)
        if strcmp(CurPaths{p}, MainPath) %Note that matlab does not save the final slash.
            HavePath = 1;
            break
        end
    end
    if HavePath == 0 %Matlab doesn't have path, so must add it
        addpath(genpath(MainPath));
    end
end

%--------------------------------------------------------------------------
%Handle various inputs from matlab or OS command lines
P = inputParser;
addOptional(P, 'InputFileT', '', @(x) ischar(x) || iscell(x) || isempty(x));
addParameter(P, 'InputFile', '', @(x) ischar(x) || iscell(x) || isempty(x));
addParameter(P, 'Chain', 'H', @(x) ismember({upper(x)}, {'H', 'L', 'HL'}));
addParameter(P, 'CheckSeqDir', 'y', @(x) ischar(x) && ismember(lower(x), {'y', 'n'}));
addParameter(P, 'Ddirection', 'all', @(x) ischar(x) && ismember(lower(x), {'all', 'fwd', 'rev', ''}));
addParameter(P, 'DevPerc', 5, @(x) isnumeric(x) && (x>=0) && (x<=100));
addParameter(P, 'Delimiter', '', @(x) ischar(x) && ismember(x, {';', ',', '\t', ''}));
addParameter(P, 'FileType', '', @ischar); %Will make input reader determine file type
addParameter(P, 'NumProc', 'max', @(x) ischar(x) || isnumeric(x));
addParameter(P, 'SettingFile', '', @ischar);
addParameter(P, 'Species', '', @(x) ischar(x) && ismember(lower(x), {'human', 'mouse', 'macaque'}));
addParameter(P, 'Strain', 'all', @ischar);
addParameter(P, 'StatusHandle', [], @(x) ishandle(x) || isempty(x) || strcmpi(class(x), 'matlab.ui.control.UIControl'));
addParameter(P, 'Vfunction', 'all', @(x) ischar(x) && min(ismember(regexpi(lower(x), ',', 'split'), {'all', 'f', 'p', 'orf', ''}))==1);
addParameter(P, 'SeqRange', [1 Inf], @(x) isnumeric(x) || ischar(x));
addParameter(P, 'Resume', 'n', @(x) ischar(x) && ismember(lower(x), {'y', 'n'})); %Resumes from known raw file
addParameter(P, 'ResumeFrom', '', @(x) ischar(x)); %File directory storing the *Raw.csv file(s).
addParameter(P, 'OutputFile', [], @(x) ischar(x) || iscell(x) || isempty(x));
addParameter(P, 'BatchSize', 1000, @(x) isnumeric(x) && x >= 1);

varargin = cleanCommandLineInput(varargin{:});
varargout = cell(1, nargout);

%Get and show the version number
if ~isempty(varargin) && ischar(varargin{1}) && ismember(varargin{1}, {'getversion', 'version'})
    if nargout == 0
        fprintf('BRILIA VERSION %s\n\n', Version);
    else
        varargout{1} = Version;
    end
    return;
end

%Determine if using a BRILIA subfunction
if ~isempty(varargin) && ~isempty(varargin{1}) && ischar(varargin{1}) && ... %Is a word
    isempty(regexpi(varargin{1}, '[^\w\d\_]')) && ... %No odd symbol
    varargin{1}(1) == lower(varargin{1}(1)) %lowerCamel    
    try
        FH = str2func(varargin{1});
        FH(varargin{2:end});
        return;
    catch ME
        rethrow(ME);
    end
end

[Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInput(P, varargin{:});
if ReturnThis
   Ps = rmfield(Ps, 'InputFileT');
   varargout = {Ps, Pu, ExpPs, ExpPu};
   return;
end

if ~isempty(Ps.InputFileT) && isempty(Ps.InputFile)
    Ps.InputFile = Ps.InputFileT;
end

%Override defaults with what is in the SettingFile
if ~isempty(Ps.SettingFile)
    Ps = readSettingFile(Ps.SettingFile, Ps);
end

%Evaluate all those that are numerical
FieldNames = fieldnames(Ps);
for j = 1:length(FieldNames)
    if ~isempty(Ps.(FieldNames{j}))         
        try 
            Ps.(FieldNames{j}) = eval(Ps.(FieldNames{j}));
        catch
        end
    end
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

%--------------------------------------------------------------------------
%Display credits, etc.
fprintf('\n');
fprintf('  Running BRILIA v%s\n', Version);
fprintf('  Developed at BHSAI\n');
fprintf('  Website: bhsai.org\n');
fprintf('\n');

%--------------------------------------------------------------------------
%Check the input and output files

%Get the full file names of input sequences
if isempty(InputFile) %Ask user to choose
    [InputFileNames, InputFilePath] = openFileDialog('*.fa*;*.*sv', 'Select the input sequence files', 'multiselect', 'on');
    if isnumeric(InputFileNames) || isempty(InputFileNames)
        return;
    end
    if ischar(InputFileNames)
        InputFileNames = {InputFileNames};
    end
    InputFile = cell(length(InputFileNames), 1);
    for f = 1:length(InputFileNames)
        InputFile{f} = [InputFilePath, InputFileNames{f}];
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
   error('%s: No valid input files were provided.', mfilename);
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
if nargin <= 1 %User did not specify species and chain, so ask user.
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
    while 1
        Selection = input('Select option: ', 's');
        if isempty(Selection)
            Selection = '1';
        end
        try 
            Selection = round(eval(Selection));
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
DB = getGeneDatabase(Species);
DB = filterGeneDatabase(DB, 'Strain', Strain, 'Ddirection', Ddirection, 'Vfunction', Vfunction);

%Set the number of processors to use
showStatus('Setting up parallel computing ...', StatusHandle);
setParallelProc(NumProc);

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

            %Check input sequence for bad characters
            [VDJdata, BadIdx] = fixInputSeq(VDJdata, VDJheader);
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
            VDJdata = seedCDR3position(VDJdata, VDJheader, DB, 'V', 80, 2, CheckSeqDir);
            VDJdata = seedCDR3position(VDJdata, VDJheader, DB, 'J', 3, 14, 'n');
            VDJdata = seedCDR3position(VDJdata, VDJheader, DB, 'Vk, Vl', 80, 2, CheckSeqDir);
            VDJdata = seedCDR3position(VDJdata, VDJheader, DB, 'Jk, Jl', 3, 14, 'n');

            %Search for initial VDJ alignment matches
            showStatus('Finding initial-guess V(D)J annotations ...', StatusHandle)
            [VDJdata, BadIdx] = findVDJmatch(VDJdata, VDJheader, DB, 'Update', 'Y');
            if max(BadIdx) ~= 0
                saveSeqData([TempDir ErrFileName], VDJdata(BadIdx, :), VDJheader, 'append');
                VDJdata(BadIdx, :) = [];
            end
            %Search for initial VJ alignment matches
            [VDJdata, BadIdx] = findVJmatch(VDJdata, VDJheader, DB, 'Update', 'Y');
            if max(BadIdx) ~= 0
                saveSeqData([TempDir ErrFileName], VDJdata(BadIdx, :), VDJheader, 'append');
                VDJdata(BadIdx, :) = [];
            end

            %Send all Non-functional to error
            [H, L, Chain] = getAllHeaderVar(VDJheader);
            FunctLoc = [H.FunctLoc L.FunctLoc];
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

            %Save remaining sequences to temp dir
            if SegmentMode %If sequence file is segmented by CDR3 length
                CDR3Loc = [H.CDR3Loc(1) L.CDR3Loc(1)];
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
    for t = 1:length(FileList)
        %Check if a Final.csv file already exists to skip
        TempRawFileName = FileList(t).name;
        TempOutputFileName = strrep(TempRawFileName, 'Raw.csv', 'Final.csv');
        if exist(TempOutputFileName, 'file')
            continue;
        end
        
        %Reload the sequence file with same-length CDR3s
        showStatus(sprintf('Processing %s ...', TempRawFileName), StatusHandle);
        [VDJdata, VDJheader] = openSeqData([TempDir TempRawFileName]);

        %Fix insertion/deletion in V framework
        showStatus('Fixing indels in V genes ...', StatusHandle)
        VDJdata = fixGeneIndel(VDJdata, VDJheader, DB);

        %Remove pseudogenes from degenerate annotations containing functional ones.
        showStatus('Accepting F genes over ORF/P ...', StatusHandle)
        VDJdata = fixDegenVDJ(VDJdata, VDJheader, DB);

        %Insure that V and J segments cover the CDR3 region.
        showStatus('Anchoring 104C and 118W/F ...', StatusHandle)
        VDJdata = constrainGeneVJ(VDJdata, VDJheader, DB);

        %Cluster the data based variable region and hamming dist of DevPerc%.
        showStatus('Clustering by lineage ...', StatusHandle)
        [VDJdata, BadVDJdataT] = clusterGene(VDJdata, VDJheader, DevPerc);
        if size(BadVDJdataT, 1) ~= 0
            saveSeqData([TempDir ErrFileName], BadVDJdataT, VDJheader, 'append');
            clear BadVDJdataT;
        end
        if size(VDJdata, 1) == 0; continue; end
        
        %Renumbering groups
        H = getHeavyHeaderVar(VDJheader);
        GrpNums = cell2mat(VDJdata(:, H.GrpNumLoc)) + GrpNumStart;
        VDJdata(:, H.GrpNumLoc) = num2cell(GrpNums);
        GrpNumStart = max(GrpNums); 
        
        %Set all groups to have same annotation and VMDNJ lengths.
        showStatus('Correcting annotations by lineage ...', StatusHandle)
        VDJdata = conformGeneGroup(VDJdata, VDJheader, DB);

        %Get better D match based on location of consensus V J mismatches.
        showStatus('Refining D annotations ...', StatusHandle)
        VDJdata = findBetterD(VDJdata, VDJheader, DB);

        %Trim V, D, J edges and extract better N regions
        showStatus('Trimming N regions ...', StatusHandle)
        VDJdata = trimGeneEdge(VDJdata, VDJheader, DB);

        %Fix obviously incorrect trees.
        showStatus('Rerooting lineage trees ...', StatusHandle)
        VDJdata = fixTree(VDJdata, VDJheader);

        %Finalize VDJdata details and CDR 1, 2, 3 info
        VDJdata = padtrimSeqGroup(VDJdata, VDJheader, 'grpnum', 'trim', 'Seq'); %will only remove "x" before and after Seq if they all have it. 
        VDJdata = findCDR1(VDJdata, VDJheader, DB);
        VDJdata = findCDR2(VDJdata, VDJheader, DB);
        VDJdata = findCDR3(VDJdata, VDJheader, DB, 'IMGT'); %removes the 104C and 118W from CDR3, and adjusts the CDR3 length to true IMGT length

        %Move non-function sequences to Err file too
        [H, L, ~] = getAllHeaderVar(VDJheader);
        FunctLoc = [H.FunctLoc(:); L.FunctLoc(:)];
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
        [VDJdata, VDJheader2] = buildVDJalignment(VDJdata, VDJheader, DB); %Adds the alignment information
        saveSeqData([TempDir TempOutputFileName], VDJdata, VDJheader2, 'append');
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

%Tries to evalues numbers, and remove leading dashes in field names (eg,
%EX:  '-saveas' => 'saveas'
%EX:  '[1,3]'   =>  [1 3]
function varargin = cleanCommandLineInput(varargin)
%Basic trial-error interpretation of number inputs
for j = 1:length(varargin)
    if ischar(varargin{j}) && ~isempty(varargin{j})
        InputStr = varargin{j};
        if isempty(regexpi(InputStr, '[^0-9\,\.\:\[\]\(\)]', 'once'))
            BrackLoc = regexpi(InputStr, '\[\]\(\)');
            if ~(InputStr(1) == '.') %In case users inputs a '.3' prefix as a string, skip.
                InputStr(BrackLoc) = [];
                try
                    varargin{j} = eval(['[' InputStr ']']);
                    continue;
                catch
                end
            end
        end
        
        if InputStr(1) == '-'
            NonDashLoc = 1;
            while InputStr(NonDashLoc) == '-'
                NonDashLoc = NonDashLoc + 1;
            end
            varargin{j} = InputStr(NonDashLoc:end);
            continue;
        end
    end
end
