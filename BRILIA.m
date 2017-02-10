%BRILIA is the main code that will run the annotation of immunosequencing
%data. BRILIA stands for B-cell Repertoire Inductive Lineage and
%Immunosequence Annotator, where "immunosequence" is simply referring to a
%DNA sequences of a B cell receptor. BRILIA works better if the sequences
%are obtained from a single repertoire of B cells from the same host.
%
%  BRILIA(FullFileNames,SettingName1,SettingValue1,...)
%
%  [RunTime] = BRILIA(FullFileNames,...)
%
%  [RunTime, SeqCount] = BRILIA(FullFileNames,...)
%
%  INPUT
%    FullFileNames: char or cell name(s) of the sequence file(s) to process
%    Param, Value pairs are for specific setting as follows
%         Setting Name    Valid Options           This sets:
%         ------------    ----------------------  -------------------------
%         SettingFile    [SettingsFile.txt]      All settings specified by
%                                                   a txt file
%
%         Species         'human' 'mouse' etc     VDJ database by species
%         Strain          'all' 'C57BL' etc       VDJ database by strain
%         Ddirection      'all' 'fwd' 'inv'       Allowed D gene direction
%         Vfunction       'all' 'f' 'p' 'orf'     Allowed V gene functions
%         DevPerc         0 <= N <= 100           Clustering cutoff dist by
%                                                   of seq length
%         FileType        'fasta', 'fastaq',      Input file type
%                         'excel', 'delimited'    
%         Delimiter       ';' ',' '\t' ''         Delimiter type (empty for
%                                                   excel file formats)
%         CheckSeqDir     'y' 'n'                 Whether or not to check 
%                                                   for complement seq
%         NumProc         'max' or N              Number of processors to
%                                                   use. Default 'max'.
%         StatusHandle    text handle obj         Handle from GUI used to 
%                                                   update status text.
%
%  OUTPUT
%    RunTime: the time, in seconds, for completing the annotation
%    SeqCount: the number of sequence that was processed
%
%  NOTE
%    The output annotation file will be saved automatically in the folder
%    where the sequence input files are, in a new folder called BRILIA.
%
%  VERSION 2.0.7
%    Version numbering is X.Y.Z, where
%      X increments for changes that affect version compatibility
%      Y increments for added/removed features, algorithm changes
%      Z increments for bug fixes and code cleaning updates
%      
%  Written by Donald Lee, dlee@bhsai.org
%  Last updated on Jan 2017

function varargout = BRILIA(varargin)
Version = '2.0.7'; %Version Number

%--------------------------------------------------------------------------
%Make sure BRILIA paths are added correctly

CurPaths = regexp(path,';','split')';
MainPath = mfilename('fullpath');
SlashLoc = regexp(MainPath,'\\|\/');
MainPath = MainPath(1:SlashLoc(end)-1);

HavePath = 0;
for p = 1:length(CurPaths)
    if strcmp(CurPaths{p},MainPath(1:end)) %Note that matlab does not save the final slash.
        HavePath = 1;
        break
    end
end
if HavePath == 0 %Matlab doesn't have path, so must add it
    disp('Adding BRILIA paths to MatLab');
    addpath(genpath(MainPath));
end

%--------------------------------------------------------------------------
%Parse the inputs

%The special case check for extracting input fields
JustGettingInput = 0;
if nargin == 1 && ischar(varargin{1})
    if strcmpi(varargin{1},'getinput')
        JustGettingInput = 1;
        varargin = {}; %Want to get defaults.
    end
end

%Check and see if there is a structure as an input. Bypass parsing then.
HasStructInput = 0;
for j = 1:length(varargin)
    if isstruct(varargin{j})
        HasStructInput = 1;
        P = varargin{j};
        break
    end
end
if HasStructInput == 0
    P = inputParser;
    addOptional(P,'FullFileNames',[],@(x) ischar(x) || iscell(x) || isempty(x));
    addParameter(P,'Species','',@(x) ischar(x) && ismember(lower(x),{'human','mouse',''})); %Note, anything empty for DB filter will prompt user input.
    addParameter(P,'Strain','',@ischar);
    addParameter(P,'Ddirection','all',@(x) ischar(x) && ismember(lower(x),{'all','fwd','inv',''}));
    addParameter(P,'Vfunction','all',@(x) ischar(x) && min(ismember(regexpi(lower(x),',','split'),{'all','f','p','orf',''}))==1);
    addParameter(P,'DevPerc',5,@(x) isnumeric(x) && (x>=0) && (x<=100)); %For clustering purposes. Set empty to force user to input it later.
    addParameter(P,'FileType','',@ischar); %Will make input reader determine file type
    addParameter(P,'Delimiter',';',@(x) ischar(x) && ismember(x,{';' ',' '\t' ''}));
    addParameter(P,'CheckSeqDir','n',@ischar);
    addParameter(P,'SettingFile','',@ischar);
    addParameter(P,'NumProc',[],@(x) ischar(x) || isnumeric(x));
    addParameter(P,'StatusHandle',[],@(x) ischar(x) || isnumeric(x));
    parse(P,varargin{:});
    P = P.Results; %For readability and for readSettingFile, remove the middle Results field.
    
    %Return inputs and stop here
    if JustGettingInput == 1
        varargout{1} = P;
        return
    end
    
    %Override defaults with what is in the SettingFile
    if ~isempty(P.SettingFile)
        P = readSettingFile(P.SettingFile,P); %You input P here to have the txt file override values in P.
    end
end

%--------------------------------------------------------------------------
%Prepare file names and VDJ databases

%Load databases and filter reference genese according to specifications
[Vmap, Dmap, Jmap] = getCurrentDatabase('change',P.Species);
[Vmap, Dmap, Jmap] = filterRefGene(Vmap,Dmap,Jmap,'Strain',P.Strain,'Ddirection',P.Ddirection,'Vfunction',P.Vfunction,'KeepThis','yes');

%Get the file names
if isempty(P.FullFileNames)
    [FileNames,FilePath] = uigetfile('*.fa*;*.xls*;*.csv;*.tsv','Select the input sequence files','multiselect','on');
    if ischar(FileNames)
        FileNames = {FileNames};
    end
    P.FullFileNames = cell(length(FileNames),1);
    for f = 1:length(FileNames)
        P.FullFileNames{f} = [FilePath, FileNames{f}];
    end
elseif ischar(P.FullFileNames)
    P.FullFileNames = {P.FullFileNames};
end

%--------------------------------------------------------------------------
%Setup parallel computing now

%Determine current status, desired status, and adjust accordingly.
showStatus('Setting up parallel computing',P.StatusHandle);
if ~isempty(P.NumProc) %Changes this ONLY if user sets it.
    %Ensuring the user input for NumProc is valid and in numbers
    if ischar(P.NumProc) %If user input 'max', set NumProc to max.
        if strcmpi(P.NumProc,'max')
            P.NumProc = feature('numCores');
        else
            P.NumProc = 1;
        end
    else %If user input a number, ensure it's a reasonable value
        P.NumProc = round(P.NumProc); %Must be integer
        if P.NumProc < 1; P.NumProc = 1; end %Can't be less than 1
        if P.NumProc > feature('numCores'); P.NumProc = feature('numCores'); end %Can't be more than Max
    end
    
    %Decide what to do about processing
    ps = parallel.Settings;
    ps.Pool.AutoCreate = false; %Ensure that parfor is not automatically run.
    PoolName = gcp('nocreate');
    if isempty(PoolName) %No current pool, so open if needed
        if P.NumProc > 1
            parpool(P.NumProc);
        end
    else %Has a pool. check if NumProc differs from NumWorkers
        if PoolName.NumWorkers ~= P.NumProc %Need to reconfigure
            delete(gcp('nocreate'));
            if P.NumProc > 1 %If parallel processing desired, open it.
                parpool(P.NumProc);
            end
        end
    end
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%For debugging only. Save current variables and then run line-by-line
DebugModeOn = 0; %Turn on(1) or off(0) debug mode in checkVDJdata
% save('temp.mat')
% return
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%==========================================================================
%BRILIA processing begins here (GUI option cannot do multiple files yet).

RunTime = zeros(length(P.FullFileNames),1); %How long it takes per file
SeqCount = zeros(length(P.FullFileNames),1); %How many sequences per file
SaveFileNames = cell(length(P.FullFileNames),1); %The output file locations
for f = 1:length(P.FullFileNames)
    tic

    %Open file and extract nucleotide information, or directly use input NTseq
    [VDJdata,VDJheader,FileName,FilePath] = convertInput2VDJdata(P.FullFileNames{f},'FileType',P.FileType,'Delimiter',P.Delimiter);
    BadVDJdata = {}; %For storing unprocessed sequences

    %Check input sequence for bad characters
    showStatus('Removing ambiguous nucletides and odd sequences.',P.StatusHandle);
    [VDJdata,BadIdx] = fixInputSeq(VDJdata,VDJheader);
        BadVDJdata = [BadVDJdata; VDJdata(BadIdx,:)];
        VDJdata(BadIdx,:) = [];
        if isempty(VDJdata); continue; end %didn't open file right

    %Find potential CDR3 start and end locations using V and J gene seed
    %alignment. Do this here, and not when doing VDJ alignment, because users
    %might have complement sequences which must be flipped.
    showStatus('Determining sequence direction and CDR3 areas.',P.StatusHandle)
    VDJdata = seedCDR3position(VDJdata,VDJheader,Vmap,'V',15,2,P.CheckSeqDir);
    VDJdata = seedCDR3position(VDJdata,VDJheader,Jmap,'J',3,14,'n');

    %Search for initial VDJ alignment matches (no try = 4x faster)
    showStatus('Finding initial-guess VDJ annotations.',P.StatusHandle)
    [VDJdata,BadIdx] = findVDJmatch(VDJdata,VDJheader,Vmap,Dmap,Jmap,'update'); %Need to implement J's are not overrride from above
        BadVDJdata = [BadVDJdata; VDJdata(BadIdx,:)];
        VDJdata(BadIdx,:) = [];
        if isempty(VDJdata); continue; end %didn't open file right

    %Fix insertion/deletion in V framework
    showStatus('Fixing indels within V segment.',P.StatusHandle)
    VDJdata = fixGeneIndel(VDJdata,VDJheader,Vmap,Dmap,Jmap);
    checkVDJdata(VDJdata,VDJheader,'fixGeneIndel',DebugModeOn);

    %Remove pseudogenes from degenerate annotations containing functional ones.
    showStatus('Removing pseudo and ORF genes if functional genes are available.',P.StatusHandle)
    VDJdata = fixDegenVDJ(VDJdata,VDJheader,Vmap,Dmap,Jmap);
    checkVDJdata(VDJdata,VDJheader,'fixDegenVDJ',DebugModeOn);

    %Insure that V and J segments cover the CDR3 region.
    showStatus('Checking if V and J segments includes 104C and 118W.',P.StatusHandle)
    VDJdata = constrainGeneVJ(VDJdata,VDJheader,Vmap,Dmap,Jmap);
    checkVDJdata(VDJdata,VDJheader,'constrainGeneVJ',DebugModeOn);

    %Pad sequences CDR3 length also have same Seq Length (required for cluster)
    [VDJdata, BadVDJdataT] = padtrimSeqGroup(VDJdata,VDJheader,'cdr3length','max','Seq');
        BadVDJdata = [BadVDJdata; BadVDJdataT];
        clear BadVDJdataT;
        if isempty(VDJdata); continue; end %didn't open file right
    checkVDJdata(VDJdata,VDJheader,'padtrimSeqGroup',DebugModeOn);

    %Remove duplicate VDJdata entries that can cause error in tree clustering
    VDJdata = removeDupSeq(VDJdata,VDJheader);
    checkVDJdata(VDJdata,VDJheader,'removeDupSeq',DebugModeOn);

    %Cluster the data based variable region and hamming dist of P.DevPerc%.
    showStatus('Performing lineage tree clustering.',P.StatusHandle)
    VDJdata = clusterGene(VDJdata,VDJheader,P.DevPerc);
    checkVDJdata(VDJdata,VDJheader,'clusterGene',DebugModeOn);

    %Set all groups to have same annotation and VMDNJ lengths.
    showStatus('Conforming VDJ annotations within clusters.',P.StatusHandle)
    VDJdata = conformGeneGroup(VDJdata,VDJheader,Vmap,Dmap,Jmap);
    checkVDJdata(VDJdata,VDJheader,'conformGeneGroup',DebugModeOn);

    %Get better D match based on location of consensus V J mismatches.
    showStatus('Refining D annotations within clusters',P.StatusHandle)
    VDJdata = findBetterD(VDJdata,VDJheader,Vmap,Dmap,Jmap);
    checkVDJdata(VDJdata,VDJheader,'findBetterD',DebugModeOn);

    %Trim V, D, J edges and extract better N regions
    showStatus('Refining N regions within clusters by trimming VDJ',P.StatusHandle)
    VDJdata = trimGeneEdge(VDJdata,VDJheader);
    checkVDJdata(VDJdata,VDJheader,'trimGeneEdge',DebugModeOn);

    %Fix obviously incorrect trees.
    showStatus('Fixing obvious errors in lineage trees.',P.StatusHandle)
    VDJdata = fixTree(VDJdata,VDJheader);
    checkVDJdata(VDJdata,VDJheader,'fixTree',DebugModeOn);

    %Finalize VDJdata details
    VDJdata = padtrimSeqGroup(VDJdata,VDJheader,'grpnum','trim','RefSeq'); %will only remove "x" before and after RefSeq. In case input seq is longer than ref seq.
    VDJdata = padtrimSeqGroup(VDJdata,VDJheader,'grpnum','trim','Seq'); %will only remove "x" before and after Seq if they all have it. 
    VDJdata = buildRefSeq(VDJdata,VDJheader,'germline','first'); %must do first seq of all cluster only
    VDJdata = updateVDJdata(VDJdata,VDJheader,Vmap,Dmap,Jmap);

    %==========================================================================
    %Save the final annotation data

    %Make a new save folder and save name
    SaveDelimiter = ';';
    SlashType = FilePath(end);
    SavePath = [FilePath 'BRILIA' SlashType];
    if ~exist(SavePath,'dir')
        mkdir(SavePath);
    end
    DotLoc = find(FileName == '.');
    if isempty(DotLoc)
        DotLoc = length(FileName);
    end

    %Save the annotated file
    SaveFullName = sprintf('%s%s.BRILIAv%s.%s',SavePath,FileName(1:DotLoc(end)-1),Version,'csv');
    saveSeqData(SaveFullName,VDJdata,VDJheader,'Delimiter',SaveDelimiter);
    SaveFileNames{f} = SaveFullName;

    %Save any unprocessed sequence
    if ~isempty(BadVDJdata)
        BadSaveFullName = sprintf('%s%s.BRILIAv%s.Err.%s',SavePath,FileName(1:DotLoc(end)-1),Version,'csv');
        saveSeqData(BadSaveFullName,BadVDJdata,VDJheader,'Delimiter',SaveDelimiter);
    end

%     %Save the settings file
%     SaveFullName = sprintf('%s%s.SettingFile.%s',SavePath,FileName(1:DotLoc(end)-1),'txt');
%     Psave = P; %Create new structure to prevent overriding P.FullFileNames
%     Psave.FullFileNames = P.FullFileNames{f}; %Select only fth name since makeSettingFile is made per file
%     makeSettingFile(SaveFullName,Psave);

    RunTime(f) = toc;
    SeqCount(f) = size(VDJdata,1) + size(BadVDJdata,1);

    %Clear big data cells before running next file
    clear VDJdata BadVDJdata
end

%Return output if needed
if nargout >= 1
    varargout{1} = RunTime;
    if nargout >= 2
        varargout{2} = SeqCount;
        if nargout >= 3
            varargout{3} = SaveFileNames;
        end
    end
end

%Simple function for showing status of BRILIA if there is a GUI or not.
function showStatus(Msg,TextHandle)
if isempty(TextHandle)
    disp(Msg);
else
    set(TextHandle,'String',Msg);
    pause(0.01); %Required short pause to update gui. Otherwise, will appear frozen.
end