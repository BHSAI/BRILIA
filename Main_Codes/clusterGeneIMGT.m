%clusterGene will locate the 104C and 118W positio of each sequence, and
%attempt to regroup based on similarity matching within the middle CDR3
%region. 
 
function  VDJdata = clusterGeneIMGT(VDJdata,VDJheader,varargin)
DevPerc = 0.03; %Fraction of NT's allowed to deviate from the rest.

FileName = '';
if isempty(VDJdata)
    [VDJdata,VDJheader,FileName,FilePath] = openSeqData;
end

%Extract the VDJ database
if length(varargin) == 3
    Vmap = varargin{1};
    Dmap = varargin{2};
    Jmap = varargin{3};
else
    [Vmap, Dmap, Jmap] = getCurrentDatabase('change');
end
H = getHeaderVar(VDJheader);

%Performing quality check
DelThese = zeros(size(VDJdata,1),1,'logical');
VMDNJ = cell2mat(VDJdata(:,H.LengthLoc));
[DelRow,~] = find(isnan(VMDNJ));
DelRow = unique(DelRow);
DelThese(DelRow) = 1;

for j = 1:size(VDJdata,1)
    FunctStatus = VDJdata{j,H.FunctLoc};
    if ischar(FunctStatus)
        if ~strcmpi(FunctStatus,'Y');
            DelThese(j) = 1;
        end
    end
end

VMDNJsum = sum(VMDNJ,2);
VMDNJsumMode = mode(VMDNJsum);
DelThese(VMDNJsum ~= VMDNJsumMode(1))=1;
VDJdata(DelThese,:) = [];

%1) Break down clusters in to smaller size based on 104C location. We avoid
%using CDR3 for those non-productive CDR3s.

%Figure out sequences with same 104C locations
ClustMap = zeros(size(VDJdata,1),6);%Map is DataRowNum | SeqNum | LeftLength | Vfamily | Dfamily | Jfamily
FamCell = VDJdata(:,H.FamNumLoc);
for j = 1:size(VDJdata,1)
    for k = 1:3
        FamCell{j,k} = FamCell{j,k}(1);
    end
end
ClustMap(:,4:6) = cell2mat(FamCell);
ClustMap(:,1) = 1:size(ClustMap,1);
ClustMap(:,2) = cell2mat(VDJdata(:,H.SeqNumLoc));
ClustMap(:,3) = cell2mat(VDJdata(:,H.CDR3Loc(2)));

%Reorganize VDJdata and clusters
ClustMap = sortrows(ClustMap,[3 4 5]);
VDJdata = VDJdata(ClustMap(:,1),:);
ClustMap(:,1) = 1:size(VDJdata,1);

%Perform crude clustering on V's to reduce clustering overload
MaxDev = ceil(length(VDJdata{1,H.SeqLoc})*DevPerc);
[~,~,UnqComboIdx] = unique(ClustMap(:,3:5),'rows');
GrpNumCt = 0;
for j = 1:max(UnqComboIdx)
    ClustIdx = ClustMap(UnqComboIdx == j,1);
    Tdata = VDJdata(ClustIdx,:);

    %Cluster by lineage
    [AncMapCell,Tdata] = buildTreeLink(Tdata,VDJheader,MaxDev);

    %Adjust the group numbers
    GrpNum2 = cell2mat(Tdata(:,H.GrpNumLoc))+GrpNumCt;
    GrpNumCt = max(GrpNum2); %Update grp num starting point
    Tdata(:,H.GrpNumLoc) = num2cell(GrpNum2);
    
    VDJdata(ClustIdx,:) = Tdata;
end

if ~isempty(FileName)
    %Save the files
    DotLoc = find(FileName == '.');
    DotLoc = DotLoc(end);
    SaveName = FileName(1:DotLoc-1);
    %Before saving to xlsx, convert columns with matrix values into char
    for q = 1:size(VDJdata,1)
        for w = 1:3
            VDJdata{q,H.FamNumLoc(w)} = mat2str(VDJdata{q,H.FamNumLoc(w)});
        end
    end
    if ispc
        xlswrite([FilePath SaveName num2str(DevPerc*100) '%.xlsx'],[VDJheader; VDJdata]);
    else
        writeDlmFile([VDJheader;VDJdata],[FilePath SaveName num2str(DevPerc*100) '.csv'],'\t');
    end
end
