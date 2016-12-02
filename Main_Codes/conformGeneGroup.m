%conformGeneGroup will force all identical groups to have the same Gene
%family and VMDNJ segment lengths

function VDJdata = conformGeneGroup(VDJdata,NewHeader,varargin)
%Extract the VDJ database
if length(varargin) == 3
    Vmap = varargin{1};
    Dmap = varargin{2};
    Jmap = varargin{3};
else
    [Vmap, Dmap, Jmap] = getCurrentDatabase;
end

getHeaderVar;

%Perform 1st group therapy: Conforming VMDNJ lengths to group
GrpNum = cell2mat(VDJdata(:,GrpNumLoc));
GrpNumUnq = unique(GrpNum);
IdxMap = 1:size(GrpNum,1);
for y = 1:length(GrpNumUnq)
    %Determine if group correction is possible
    GrpLoc = (GrpNumUnq(y) == GrpNum);
    IdxLoc = IdxMap(GrpLoc);
    if sum(GrpLoc) <= 1; continue; end %Nothing to correct, single entry.

    %Standardize the other fields as well (1st one is the root)
    Tdata = VDJdata(IdxLoc,:);    
    Tdata(:,[LengthLoc FamLoc FamNumLoc DelLoc]) = repmat(Tdata(1,[LengthLoc FamLoc FamNumLoc DelLoc]),size(Tdata,1),1);

    VDJdata(IdxLoc,:) = Tdata;  
end

VDJdata = buildVDJalignment(VDJdata,NewHeader,Vmap,Dmap,Jmap);
VDJdata = makeClassifier(VDJdata,NewHeader);
VDJdata = appendMutCt(VDJdata,NewHeader); %SHM info on the VMDNJ segments
