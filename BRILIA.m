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
%         SettingsFile    [SettingsFile.txt]      All settings specified by
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
%
%  OUTPUT
%    RunTime: the time, in seconds, for completing the annotation
%    SeqCount: the number of sequence that was processed
%
%  NOTE
%    The output annotation file will be saved automatically in the folder
%    where the sequence input files are, in a new folder called BRILIA.
%
%  VERSION 2.0.2
%    Version numbering is X.Y.Z, where
%      X increments for substantial changes that alter output results
%      Y increments for added/removed features
%      Z increments for bug fixes and code cleaning updates
%      
%  Written by Donald Lee, dlee@bhsai.org
%  Last updated on Jan 2017

function varargout = BRILIA(varargin)
Version = '2.0.2'; %Version Number

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

%If settings file is given, use it
HaveSettingFile = 0;
if ~isempty(varargin)
    SettingLoc = findCell('SettingsFile',varargin,'MatchCase','any');
    if min(SettingLoc) > 0
        SettingFileName = varargin{SettingLoc+1};
        P = readSettingFile(SettingFileName);
        HaveSettingFile = 1;
    end
end

%Otherwise, use default + user inputs
if HaveSettingFile == 0
    P = inputParser;
    addOptional(P,'FullFileNames',[],@(x) ischar(x) || iscell(x) || isempty(x));
    addParameter(P,'Species','',@(x) ischar(x) && ismember(lower(x),{'human','mouse'})); %Note, anything empty for DB filter will prompt user input.
    addParameter(P,'Strain','',@ischar);
    addParameter(P,'Ddirection','all',@(x) ischar(x) && ismember(lower(x),{'all','fwd','inv'}));
    addParameter(P,'Vfunction','all',@(x) ischar(x) && min(ismember(regexpi(lower(x),',','split'),{'all','f','p','orf'}))==1);
    addParameter(P,'DevPerc',3,@(x) isnumeric(x) && (x>=0) && (x<=100)); %For clustering purposes. Set empty to force user to input it later.
    addParameter(P,'FileType','',@ischar); %Will make input reader determine file type
    addParameter(P,'Delimiter',';',@(x) ischar(x) && ismember(x,{';' ',' '\t' ''}));
    addParameter(P,'CheckSeqDir','n',@ischar);
    parse(P,varargin{:});
    P = P.Results; %For readability, remove the middle Results field.
end

FullFileNames = P.FullFileNames;
DevPerc = P.DevPerc;
Vfunction = P.Vfunction;
Ddirection = P.Ddirection;
Species = P.Species;
Strain = P.Strain;
FileType = P.FileType;
Delimiter = P.Delimiter;
CheckSeqDir = P.CheckSeqDir;

%Load databases and filter reference genese according to specifications
[Vmap, Dmap, Jmap] = getCurrentDatabase('change',Species);
[Vmap, Dmap, Jmap] = filterRefGene(Vmap,Dmap,Jmap,'Strain',Strain,'Ddirection',Ddirection,'Vfunction',Vfunction,'KeepThis','yes');

%Get the file names
if isempty(FullFileNames)
    [FileNames,FilePath] = uigetfile('*.fa*;*.xls*;*.csv;*.tsv','Select the input sequence files','multiselect','on');
    if ischar(FileNames)
        FileNames = {FileNames};
    end
    FullFileNames = cell(length(FileNames),1);
    for f = 1:length(FileNames)
        FullFileNames{f} = [FilePath, FileNames{f}];
    end
elseif ischar(FullFileNames)
    FullFileNames = {FullFileNames};    
end

RunTime = zeros(length(FileNames),1);
SeqCount = zeros(length(FileNames),1);
 
for f = 1:length(FullFileNames)
    try
        tic

        %Open file and extract nucleotide information, or directly use input NTseq
        [VDJdata,NewHeader,FileName,FilePath] = convertInput2VDJdata(FullFileNames{f},'FileType',FileType,'Delimiter',Delimiter);
        BadVDJdata = {}; %For storing unprocessed sequences
        getHeaderVar;

        %==========================================================================
        %BRILIA processing begins here

        %Configure the parallel processing feature for larger datasets
        ps = parallel.Settings;
        ps.Pool.AutoCreate = false; %Ensure that parfor is not automatically run.
        PoolName = gcp('nocreate');
        if size(VDJdata,1) > 200 && isempty(PoolName) == 1
            CoreNum = feature('numCores');
            parpool(CoreNum);
        end

        %Check input sequence for bad characters
        disp('Removing ambiguous nucletides and odd sequences.')
        [VDJdata,BadIdx] = fixInputSeq(VDJdata,NewHeader);
            BadVDJdata = [BadVDJdata; VDJdata(BadIdx,:)];
            VDJdata(BadIdx,:) = [];

        %Find potential CDR3 start and end locations using V and J gene seed
        %alignment. Do this here, and not when doing VDJ alignment, because users
        %might have complement sequences which must be flipped.
        disp('Determining sequence direction and CDR3 areas.')
        VDJdata = seedCDR3position(VDJdata,NewHeader,Vmap,'V',15,2,CheckSeqDir);
        VDJdata = seedCDR3position(VDJdata,NewHeader,Jmap,'J',3,14,'n');

        %Search for initial VDJ alignment matches (no try = 4x faster)
        disp('Finding initial-guess VDJ annotations.')
        [VDJdata,BadIdx] = findVDJmatch(VDJdata,NewHeader,Vmap,Dmap,Jmap,'update'); %Need to implement J's are not overrride from above
            BadVDJdata = [BadVDJdata; VDJdata(BadIdx,:)];
            VDJdata(BadIdx,:) = [];

        %For debugging only. Save current variables so you can run each code,
        %line by line.
        save('temp.mat')
        checkVDJdata(VDJdata,NewHeader,'findVDJmatch');

        %Fix insertion/deletion in V framework
        disp('Fixing indels within V segment.')
        VDJdata = fixGeneIndel(VDJdata,NewHeader,Vmap,Dmap,Jmap);
        checkVDJdata(VDJdata,NewHeader,'fixGeneIndel');

        %Remove pseudogenes from degenerate annotations containing functional ones.
        disp('Removing pseudo and ORF genes if functional genes are available.')
        VDJdata = fixDegenVDJ(VDJdata,NewHeader,Vmap,Dmap,Jmap);
        checkVDJdata(VDJdata,NewHeader,'fixDegenVDJ');

        %Insure that V and J segments cover the CDR3 region.
        disp('Checking if V and J segments includes 104C and 118W.')
        VDJdata = constrainGeneVJ(VDJdata,NewHeader,Vmap,Dmap,Jmap);
        checkVDJdata(VDJdata,NewHeader,'constrainGeneVJ');

        %Pad sequences CDR3 length also have same Seq Length (required for cluster)
        [VDJdata, BadVDJdataT] = padtrimSeqGroup(VDJdata,NewHeader,'cdr3length','max');
            BadVDJdata = [BadVDJdata; BadVDJdataT];
            clear BadVDJdataT;
        checkVDJdata(VDJdata,NewHeader,'padtrimSeqGroup');

        %Remove duplicate VDJdata entries that can cause error in tree clustering
        VDJdata = removeDupSeq(VDJdata,NewHeader);
        checkVDJdata(VDJdata,NewHeader,'removeDupSeq');

        %Cluster the data based variable region and hamming dist of DevPerc%.
        disp('Performing lineage tree clustering.')
        VDJdata = clusterGene(VDJdata,NewHeader,DevPerc);
        checkVDJdata(VDJdata,NewHeader,'clusterGene');

        %Set all groups to have same annotation and VMDNJ lengths.
        disp('Conforming VDJ annotations within clusters.')
        VDJdata = conformGeneGroup(VDJdata,NewHeader,Vmap,Dmap,Jmap);
        checkVDJdata(VDJdata,NewHeader,'conformGeneGroup');

        %Get better D match based on location of consensus V J mismatches.
        disp('Refining D annotations within clusters')
        VDJdata = findBetterD(VDJdata,NewHeader,Vmap,Dmap,Jmap);
        checkVDJdata(VDJdata,NewHeader,'findBetterD');

        %Trim V, D, J edges and extract better N regions
        disp('Refining N regions within clusters by trimming VDJ')
        VDJdata = trimGeneEdge(VDJdata,NewHeader);
        checkVDJdata(VDJdata,NewHeader,'trimGeneEdge');

        %Fix obviously incorrect trees.
        disp('Fixing obvious errors in lineage trees.')
        VDJdata = fixTree(VDJdata,NewHeader);
        checkVDJdata(VDJdata,NewHeader,'fixTree');

        %Finalize VDJdata details
        VDJdata = padtrimSeqGroup(VDJdata,NewHeader,'grpnum','trim'); %will only remove "x" before and after sequences if they all have it. 
        VDJdata = buildRefSeq(VDJdata,NewHeader,'germline','first'); %must do first seq of all cluster only
        VDJdata = updateVDJdata(VDJdata,NewHeader,Vmap,Dmap,Jmap);

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
        saveSeqData(SaveFullName,VDJdata,NewHeader,'Delimiter',SaveDelimiter);

        %Save any unprocessed sequence
        if ~isempty(BadVDJdata)
            BadSaveFullName = sprintf('%s%s.BRILIAv%s.Err.%s',SavePath,FileName(1:DotLoc(end)-1),Version,'csv');
            saveSeqData(BadSaveFullName,BadVDJdata,NewHeader,'Delimiter',SaveDelimiter);
        end

        RunTime(f) = toc;
        SeqCount(f) = size(VDJdata,1) + size(BadVDJdata,1);

        %Clear big data cells before running next file
        clear VDJdata BadVDJdata
    catch
        %Clear big data cells before running next file
        disp(['Fatal error processing: ' FileNames{f} ' . Check Delimiter and FileType inputs']);
        clear VDJdata BadVDJdata
    end
end

%Return output if needed
if nargout >= 1
    varargout{1} = RunTime;
    if nargout >= 2
        varargout{2} = SeqCount;
    end
end
