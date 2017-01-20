%BRILIA is the core function that takes in one input file and returns one
%output file. If summoned via BRILIAbatch, then multiple input files can be
%processed.

function varargout = BRILIA(varargin)
Version = '2.0.1'; %Version Number

%--------------------------------------------------------------------------
%Parse the inputs
P = inputParser;
addOptional(P,'FullFileName',[],@(x) ischar(x) || iscell(x) || isempty(x));
addParameter(P,'Species','',@(x) ischar(x) && ismember(lower(x),{'human','mouse'})); %Note, anything empty for DB filter will prompt user input.
addParameter(P,'Strain','all',@ischar);
addParameter(P,'Ddirection','all',@(x) ischar(x) && ismember(lower(x),{'all','fwd','inv'}));
addParameter(P,'Vfunction','all',@(x) ischar(x) && min(ismember(regexpi(lower(x),',','split'),{'all','f','p','orf'}))==1);
addParameter(P,'DevPerc',3,@(x) isnumeric(x) && (x>=0) && (x<=100)); %For clustering purposes. Set empty to force user to input it later.
addParameter(P,'FileType','',@ischar);
addParameter(P,'Delimiter',';',@(x) ischar(x) && ismember(x,{';' ',' '\t'}));
addParameter(P,'CheckSeqDir','y',@ischar);
%addParameter(P,'SkipFirstMatch',0,@(x) isnumeric(x) && (x==1 || x==0));
%%Skips the first findVDJmatch step, good if using BRILIA input files and
%%just want to test later steps. Debuggin only.

parse(P,varargin{:});

FullFileName = P.Results.FullFileName;
DevPerc = P.Results.DevPerc;
Vfunction = P.Results.Vfunction;
Ddirection = P.Results.Ddirection;
Species = P.Results.Species;
Strain = P.Results.Strain;
FileType = P.Results.FileType;
Delimiter = P.Results.Delimiter;
CheckSeqDir = P.Results.CheckSeqDir;
%SkipFirstMatch = P.Results.SkipFirstMatch; %If you are uploading a .raw file, you can skip to the correction part

%Used to help with debugging. Modify as needed.
% DevNum = 14;
% FullFileName = '';
% SkipFirstMatch = 0;
% DevPerc = 5;
% Vfunction = 'all';
% Ddirection = 'all';
% Species = 'mouse';
% Strain = 'C57BL';
% FileType = '';
% Delimiter = '\t';
% CheckSeqDir = 'n';

%Obtain, filter, and reduce the VDJ germgline reference maps
[Vmap, Dmap, Jmap] = getCurrentDatabase('change',Species);
[Vmap, Dmap, Jmap] = filterRefGene(Vmap,Dmap,Jmap,'Strain',Strain,'Ddirection',Ddirection,'Vfunction',Vfunction,'KeepThis','yes');

%Open file and extract nucleotide information, or directly use input NTseq
[VDJdata,NewHeader,FileName,FilePath] = convertInput2VDJdata(FullFileName,'FileType',FileType,'Delimiter',Delimiter);
BadVDJdata = {}; %For storing unprocessed sequences
getHeaderVar;

% if SkipFirstMatch == 1 %Debugging only
%     [VDJdata,NewHeader,FileName,FilePath] = openSeqData(FullFileName);
% end

%==========================================================================
%BRILIA processing begins here
tic

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
VDJdata = seedCDR3position(VDJdata,NewHeader,Vmap,'V',9,2,CheckSeqDir);
VDJdata = seedCDR3position(VDJdata,NewHeader,Jmap,'J',3,2,'n'); %Only use V to check direction, since J is too short.

%Search for initial VDJ alignment matches
disp('Finding initial-guess VDJ annotations.')
[VDJdata,BadIdx] = findVDJmatch(VDJdata,NewHeader,Vmap,Dmap,Jmap,'update'); %Need to implement J's are not overrride from above
    BadVDJdata = [BadVDJdata; VDJdata(BadIdx,:)];
    VDJdata(BadIdx,:) = [];

%Fix insertion/deletion in V framework
disp('Fixing indels within V segment.')
VDJdata = fixGeneIndel(VDJdata,NewHeader,Vmap,Dmap,Jmap);
    
%Remove pseudogenes from degenerate annotations containing functional ones.
disp('Removing pseudo and ORF genes if functional genes are available.')
VDJdata = fixDegenVDJ(VDJdata,NewHeader,Vmap,Dmap,Jmap);

%Insure that V and J segments cover the CDR3 region.
disp('Checking V and J segments covers 104C and 118W.')
VDJdata = constrainGeneVJ(VDJdata,NewHeader,Vmap,Dmap,Jmap);

%Pad sequences CDR3 length also have same Seq Length (required for cluster)
[VDJdata, BadVDJdataT] = padtrimSeqGroup(VDJdata,NewHeader,'cdr3length','max');
    BadVDJdata = [BadVDJdata; BadVDJdataT];
    clear BadVDJdataT;

%Remove duplicate VDJdata entries that can cause error in tree clustering
VDJdata = removeDupSeq(VDJdata,NewHeader);

%Cluster the data based variable region and hamming dist of DevPerc%.
disp('Performing lineage tree clustering.')
VDJdata = clusterGene(VDJdata,NewHeader,DevPerc);
%VDJdata = fixTree(VDJdata,NewHeader);

%Set all groups to have same annotation and VMDNJ lengths.
disp('Conforming VDJ annotations within clusters.')
VDJdata = conformGeneGroup(VDJdata,NewHeader,Vmap,Dmap,Jmap);

%Get better D match based on location of consensus V J mismatches.
disp('Refining D annotations within clusters')
VDJdata = findBetterD(VDJdata,NewHeader,Vmap,Dmap,Jmap);

%Trim V, D, J edges and extract better N regions
disp('Refining N regions within clusters by trimming VDJ')
VDJdata = trimGeneEdge(VDJdata,NewHeader);

%Fix obviously incorrect trees.
disp('Fixing obvious errors in lineage trees.')
VDJdata = fixTree(VDJdata,NewHeader);

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

RunTime = toc;

%Return output, if any.
if nargout >= 1
    varargout{1} = VDJdata;    
    if nargout >= 2
        varargout{2} = NewHeader;
        if nargout >= 3
            varargout{3} = RunTime;
        end
    end
end