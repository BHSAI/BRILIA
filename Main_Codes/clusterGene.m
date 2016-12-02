%clusterGene will locate the 104C and 118W position of each sequence, and
%attempt to regroup based on similarity matching within the middle CDR3
%region. 

%  VDJdata = clusterGene(VDJdata,NewHeader,DevPerc)
%  VDJdata = clusterGene(VDJdata,NewHeader,DevPerc,Vmap,Dmap,Jmap)
%  DevPerc is a fraction from 0 to 1, indicating fr of sequence length for
%  SHMdistance. EX: DevPerc = 0.03 for 125 bp seq is 4. This sets cutoff
%  percent.
function  VDJdata = clusterGene(VDJdata,NewHeader,DevPerc,varargin)
%Extract the VDJ database
if length(varargin) == 3
    Vmap = varargin{1};
    Dmap = varargin{2};
    Jmap = varargin{3};
else
    [Vmap, Dmap, Jmap] = getCurrentDatabase;
end

getHeaderVar;

%Do coarse clustering by V, D, J family number and CDR3 length + pos (must
%be trimmed and aligned according to 118W at 3' end.
ClustMap = zeros(size(VDJdata,1),5); %[RowNum SeqNum CDR3length VfamNum JfamNum]
for j = 1:size(VDJdata,1)
    %Find the location of the conserved C of every V
    VmapNum = VDJdata{j,FamNumLoc(1)};    
    JmapNum = VDJdata{j,FamNumLoc(3)};
    CDR3Len = VDJdata{j,CDR3Loc(2)};
    
    %Determine the family number
    VFamNumIdx = regexp(Vmap{VmapNum(1),4},'\d');
    VfamNum = eval(Vmap{VmapNum(1),4}(VFamNumIdx));
    JFamNumIdx = regexp(Jmap{JmapNum(1),4},'\d');
    JfamNum = eval(Jmap{JmapNum(1),4}(JFamNumIdx));
       
    ClustMap(j,:) = [j VDJdata{j,SeqNumLoc} CDR3Len VfamNum JfamNum];
end

%Reorganize VDJdata and clusters
ClustMap = sortrows(ClustMap,[3 4 5]);
VDJdata = VDJdata(ClustMap(:,1),:);
ClustMap(:,1) = 1:size(VDJdata,1); %Renumber
[~,~,UnqComboIdx] = unique(ClustMap(:,3:5),'rows'); %Identify unique coarse clusters

%For each coarse cluster, find subclusters using Tree clustering and SHM
%distance cutoff
MaxDev = ceil(length(VDJdata{1,SeqLoc})*DevPerc); %Cutoff SHM distance
GrpNumCt = 0;
for j = 1:max(UnqComboIdx)
    ClustIdx = ClustMap(UnqComboIdx == j,1);
    Tdata = VDJdata(ClustIdx,:);
    [~,Tdata] = buildTreeLink(Tdata,NewHeader,MaxDev);
        
    %Adjust the group numbers
    GrpNum2 = cell2mat(Tdata(:,GrpNumLoc)) + GrpNumCt;
    GrpNumCt = max(GrpNum2); %Update grp num starting point
    Tdata(:,GrpNumLoc) = num2cell(GrpNum2);
    
    VDJdata(ClustIdx,:) = Tdata;
end