%matchVDJ will open up either VDJdata or ADAPdata file format, and then
%process the NT sequences to extract VDJ annotations.

function varargout = BRILIA(varargin)
DevNum = 13; %Version Number

P = inputParser;
addOptional(P,'FullName',[],@ischar);
addOptional(P,'Vmap',{},@iscell);
addOptional(P,'Dmap',{},@iscell);
addOptional(P,'Jmap',{},@iscell);
addParameter(P,'SkipFirstMatch',0,@isnumeric);
addParameter(P,'DevPerc',[],@isnumeric);
parse(P,varargin{:});

FullName = P.Results.FullName;
Vmap = P.Results.Vmap;
Dmap = P.Results.Dmap;
Jmap = P.Results.Jmap;
SkipFirstMatch = P.Results.SkipFirstMatch; %If you are uploading a .raw file, you can skip to the correction part
DevPerc = P.Results.DevPerc;

%==========================================================================
%If you are doing this in batch, and want to push in a DevPerc, do it here.
if isempty(DevPerc)
    DevPerc = input('What is the clustering deviation, in fraction? Default 0.03: '); 
    if isempty(DevPerc)
        DevPerc = 0.03; %0.05 for simulation.
    elseif DevPerc >= 1 || DevPerc < 0
        error('Deviation percentage should be 1 > DevPerc >= 0');
    end
end

%==========================================================================
%Open file and extract nucleotide information, or directly use input NTseq
if isempty(FullName) %Open file
    [SampleData, SampleHeader, FileName, FilePath] = openSeqData;
else
    [SampleData, SampleHeader, FileName, FilePath] = openSeqData(FullName);    
end

%==========================================================================
%Make a new save folder for current version
SlashType = FilePath(end);
SavePath = [FilePath 'BRILIAv' num2str(DevNum) SlashType];
if ~exist(SavePath,'dir')
    mkdir(SavePath);
end

%Make the savename prefix
DotLoc = find(FileName == '.');
if isempty(DotLoc)
    DotLoc = length(FileName);
end
SaveNamePre = FileName(1:DotLoc(end)-1);

%==========================================================================
%Determine if this should be starting from the beginning, or from
%correction step.
if SkipFirstMatch == 0
    %Obtain the MAVRIC VDJdata header info for output format
    [~, ~, HeaderData] = xlsread('Headers_BRILIA.xlsx');
    NewHeader = HeaderData(2:end,1)';
    VDJdata = cell(size(SampleData,1),length(NewHeader));
    getHeaderVar;
    
    %Determine the MAVRIC SeqNum values based on SeqNames, if it exists
    SeqNameLocT = findHeader(SampleHeader,{'SeqName','SeqNum'});
    
    %If you have both, do a quick check to make sure to fille SeqName by
    %SeqNum, or vice versa.
    SeqNameLocT(SeqNameLocT == 0) = [];
    if length(SeqNameLocT) == 2
        for j = 1:size(SampleData,1)
            Snum = SampleData{j,SeqNameLocT};
            Sname = SampleData{j,SeqNameLocT};
            if isempty(Sname) && ~isempty(Snum)
                SampleData{j,SeqNameLocT(2)} = Snum;
            end
        end
    end
    
    SeqNameLocT = SeqNameLocT(1);
    HaveName = 0;
    HaveUnqName = 1;
    if SeqNameLocT > 0
        HaveName = 1;
        for k = 1:size(SampleData,SeqNameLocT)
            if isempty(SampleData{k,SeqNameLocT}) || ~isnumeric(SampleData{k,SeqNameLocT})
                HaveUnqName = 0; %Can't be unique names if it's empty or has characters
                break
            end
        end
        if HaveUnqName == 1 %Check to make sure they really are unique numbers
            UnqNames = unique(cell2mat(SampleData(:,SeqNameLocT)));
            if length(UnqNames) ~= size(SampleData,1)
                HaveUnqName = 0;
            end
        end
    else
        HaveUnqName = 0;
    end
    
    if HaveName == 1 %Has sequence names
        if HaveUnqName == 0 %but is not unique names or has characters
            VDJdata(:,SeqNumLoc) = num2cell([1:size(SampleData,1)]');
        else %Otherwise just copy over seq name
            VDJdata(:,SeqNumLoc) = SampleData(:,SeqNameLocT);
        end
        VDJdata(:,SeqNameLoc) = SampleData(:,SeqNameLocT);
    else %Have no SeqName -> Assign number to both SeqNum and SeqName
        VDJdata(:,SeqNumLoc) = num2cell([1:size(SampleData,1)]');
        VDJdata(:,SeqNameLoc) = num2cell([1:size(SampleData,1)]');
    end
    
    %Assign a unique group number for now. It will change later
    VDJdata(:,GrpNumLoc) = num2cell([1:size(SampleData,1)]');
    
    %Transfer over sequence and template count (if avail)
    SeqLocT = findHeader(SampleHeader,{'nucleotide','Seq'});
        SeqLocT(SeqLocT==0) = []; SeqLocT = SeqLocT(1);
        VDJdata(:,SeqLoc) = SampleData(:,SeqLocT);
    TemplateLocT = findHeader(SampleHeader,{'count (templates)','TemplateCount'});
        TemplateLocT(TemplateLocT==0) = [];
        if isempty(TemplateLocT) %No template information, set all to 1
            VDJdata(:,TemplateLoc) = num2cell(ones(size(VDJdata,1),1));
        else
            VDJdata(:,TemplateLoc) = SampleData(:,TemplateLocT);
        end
else
    VDJdata = SampleData;
    NewHeader = SampleHeader;
    getHeaderVar;
end
clear SampleData

if isempty(Vmap) || isempty(Dmap) || isempty(Jmap)
    [Vmap, Dmap, Jmap] = getCurrentDatabase; %Selective active full database
    [Vmap, Dmap, Jmap] = filterRefGene(Vmap,Dmap,Jmap); %Selelect host strain
end
    
%Trim Vmap sequences since you don't need it all.
MaxLen = 1; %Max length of the sequences
for j = 1:size(VDJdata,1);
    SeqLen = length(VDJdata{j,SeqLoc});
    if SeqLen > MaxLen;
        MaxLen = SeqLen;
    end
end
VuseLen = MaxLen + 25; %Use V segments up to 25 nts more than MaxLen
for v1 = 1:size(Vmap,1)
    TempNT = Vmap{v1,1};
    if length(TempNT) > VuseLen
        Vmap{v1,1} = TempNT(end-VuseLen+1:end);
    end
end

%==========================================================================
%Find VDJ matches as an initial starting point.
if SkipFirstMatch == 0
    disp('Finding initial-guess VDJ annotatations')
    
    %Search for initial VDJ alignment matches
    VDJdata = findVDJmatch(VDJdata,NewHeader,Vmap,Dmap,Jmap);

    %----------------------------------------------------------------------
    %Save for debuging only
    SaveRaw = 0;
    if SaveRaw == 1
        %Save the raw annotations now
        SaveName = [SaveNamePre '.Raw'];

        %Before saving to xlsx, convert columns with matrix values into char
        VDJdataSave = VDJdata;
        for q = 1:size(VDJdata,1)
            for w = 1:3
                VDJdataSave{q,FamNumLoc(w)} = mat2str(VDJdataSave{q,FamNumLoc(w)});
            end
        end

        %Save to excel or csv file, depending on OS
        if ispc
            xlswrite([SavePath SaveName '.xlsx'],cat(1,NewHeader,VDJdataSave));
        else
            writeDlmFile(cat(1,NewHeader,VDJdataSave),[SavePath SaveName '.csv'],'\t');
        end

        clear VDJdataSave
    end
end

%==========================================================================
%Fix insertion/deletion in V framework
disp('Fixing Indels within V Framework Segment')
VDJdata = fixGeneIndel(VDJdata,NewHeader,Vmap,Dmap,Jmap);

%Remove duplicates that could arise from indel corrections
disp('Removing duplicates that formed from Indel corrections')
VDJdata = removeDupVDJdata(VDJdata,NewHeader);

%Fix J's if 118 F/W is outside the seq lenght. Happens if you have short
%seq, if end of seq is W or F, and if you can find another J that fills
%this spot.
disp('Fixing J alignmenets and no-CDR3 errors')
VDJdata = fixJalign(VDJdata,NewHeader,Vmap,Dmap,Jmap);

%Fixes multiple annotations where some are pseudo and other are functional.
%Gets rid of the pseudo annotaitons.
VDJdata = fixPseudoVDJ(VDJdata,NewHeader,Vmap,Dmap,Jmap);

%Fix V and J lengths and dels based on 104C and 118 F/W Location
disp('Fixating V and J segments by 104C and 118W')
VDJdata = constrainGeneVJ(VDJdata,NewHeader,Vmap,Dmap,Jmap);

%==========================================================================
%PERHAPS HERE, implement trimming and alignment...   
%==========================================================================

%Recluster the data based variable region and hamming dist of 5%.
disp('Performing Tree clustering')
VDJdata = clusterGene(VDJdata,NewHeader,DevPerc,Vmap,Dmap,Jmap);

%Set all groups to have same annotation and VMDNJ lengths
disp('Conforming VDJ annotations within clusters')
VDJdata = conformGeneGroup(VDJdata,NewHeader,Vmap,Dmap,Jmap,2);

%Label nonprod CDR3 at this point
disp('Labeling nonproductive VDJs')
VDJdata = labelNonprodVDJ(VDJdata,NewHeader); %need to do this since it's possible that after conforming, you get correct reading frame and no stop codons.

%Get better D match based on location of consensus V J mismatches.
disp('Refining D annotations within clusters')
VDJdata = findBetterD(VDJdata,NewHeader,Vmap,Dmap,Jmap);

%Trim V, D, J edges and extract better N regions
disp('Refining N regions within clusters by trimming VDJ')
VDJdata = trimGeneEdge(VDJdata,NewHeader,Vmap,Dmap,Jmap);

VDJdata = buildRefSeq(VDJdata,NewHeader,'same','germline','first');%Same length, germline substitution, on first sequence of each group
VDJdata = appendMutCt(VDJdata,NewHeader); %SHM infor on the VMDNJ segments

%Final lineage tree correction: At this point, it's obvious if you got the
%wrong root sequence. Redraw all tree using true nearest neighbor.
VDJdata = fixTree(VDJdata,NewHeader);

%==========================================================================
%Save the final annotation
%--------------------------------------------------------------------------
SaveName = [SaveNamePre '.BRILIA'];

%Before saving to xlsx, convert columns with matrix values into char
for q = 1:size(VDJdata,1)
    for w = 1:3
        VDJdata{q,FamNumLoc(w)} = mat2str(VDJdata{q,FamNumLoc(w)});
    end
end

%Save to excel or csv file, depending on OS
if ispc
    xlswrite([SavePath SaveName '.xlsx'],cat(1,NewHeader,VDJdata));
else
    writeDlmFile(cat(1,NewHeader,VDJdata),[SavePath SaveName '.csv'],'\t');
end
 
if nargout >= 1
    varargout{1} = VDJdata;    
    if nargout >= 2
        varargout{2} = NewHeader;
    end
end