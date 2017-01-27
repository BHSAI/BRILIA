%clusterGene will perform the lineage-based clustering scheme used by
%BRILIA. A crude cluster is first assembled based on save V and J "family
%gene numbers", and afterwards, will determine which cluster belongs
%together based on the SHM distance calculations. The cutoff distance has
%been specfied as a % of the sequence length (DevPerc), though SHM distance
%can return scores larger/smaller that DevPerc length of Seq depending on
%if the SHM patterns match with expected patterns.
%
%  VDJdata = clusterGene(VDJdata,NewHeader,DevPerc)
%
%  INPUT
%    DevPerc: Ranges 0 to 100, indicating % of sequence length to use as
%    the cutoff SHM distance.
%
%  See also buildTreeLink

function  VDJdata = clusterGene(VDJdata,NewHeader,DevPerc)
getHeaderVar;

%Create the clusters
ClustMap = zeros(size(VDJdata,1),4); %[RowNum CDR3length VfamNum JfamNum]
for j = 1:size(VDJdata,1)
    %Find length of the CDR3nts
    CDR3start = VDJdata{j,CDR3Loc(3)};
    if isempty(CDR3start); CDR3start = 0; end
    CDR3end = VDJdata{j,CDR3Loc(4)};
    if isempty(CDR3end); CDR3end = 0; end
    CDR3Len = CDR3end - CDR3start;
    
    %Determine the V family number
    Vname = VDJdata{j,FamLoc(1)};    
    Vnums = regexp(Vname,'[^\d]*','split');
    for k = 1:length(Vnums)
        if ~isempty(Vnums{k}); break; end
    end
    Vfamily = eval(Vnums{k});
    
    %Determine the J family number
    Jname = VDJdata{j,FamLoc(3)};
    Jnums = regexp(Jname,'[^\d]*','split');
    for k = 1:length(Jnums)
        if ~isempty(Jnums{k}); break; end
    end
    Jfamily = eval(Jnums{k});
       
    ClustMap(j,:) = [j CDR3Len Vfamily Jfamily];
end

%Resort VDJdata and determine crude clusters
ClustMap = sortrows(ClustMap,[2 3 4]);
VDJdata = VDJdata(ClustMap(:,1),:);
ClustMap(:,1) = 1:size(VDJdata,1); %Renumber
[~,~,CrudeClustIdx] = unique(ClustMap(:,2:4),'rows'); %Identify unique coarse clusters

%For each coarse cluster, find subclusters using Tree clustering and SHM
%distance cutoff
MaxDev = ceil(length(VDJdata{1,SeqLoc})*DevPerc/100); %Cutoff SHM distance
GrpNumCt = 0; %buildTreeLink will provide a group number clust from 1 to N
for j = 1:max(CrudeClustIdx)
    try
        ClustIdx = ClustMap(CrudeClustIdx == j,1);
        Tdata = VDJdata(ClustIdx,:);       
        [~,Tdata] = buildTreeLink(Tdata,NewHeader,MaxDev);

        %Adjust the group numbers
        GrpNum2 = cell2mat(Tdata(:,GrpNumLoc)) + GrpNumCt; %Get new grp nums
        GrpNumCt = max(GrpNum2); %Update grp num starting point
        Tdata(:,GrpNumLoc) = num2cell(GrpNum2);
    
        VDJdata(ClustIdx,:) = Tdata;
    catch
        WarningMsg = sprintf('Warning at %s, cluster group # %d',mfilename,j);
        disp(WarningMsg); pause
    end
end