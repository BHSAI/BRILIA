%clusterGene will perform the lineage-based clustering scheme used by
%BRILIA. A crude cluster is first assembled based on save V and J "family
%gene numbers", and afterwards, will determine which cluster belongs
%together based on the SHM distance calculations. The cutoff distance has
%been specfied as a % of the sequence length (DevPerc), though SHM distance
%can return scores larger/smaller that DevPerc length of Seq depending on
%if the SHM patterns match with expected patterns.
%
%  VDJdata = clusterGene(VDJdata, VDJheader, DevPerc)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    VDJheader: main BRILIA header cell
%    DevPerc: Ranges 0 to 100, indicating % of sequence length to use as
%      the cutoff SHM distance.
%
%  OUTPUT
%    VDJdata: modified VDJdata where sequences are clustered based on
%      lineage relationships and maximum SHM distance, computed as:
%      CutoffDist = DevPerc x SeqLength
%
%  See also buildTreeLink

function [VDJdata, BadVDJdata] = clusterGene(VDJdata, Map, DevPerc)
%Determine chain and extract key locations
switch Map.Chain
    case 'HL'
        FunctIdx = [Map.hFunct Map.lFunct];
        ExtractIdx = [Map.hCDR3(2) Map.hGeneName(1) Map.hGeneName(end) Map.lCDR3(2) Map.lGeneName(1) Map.lGeneName(end)];
        SortOrder = [1 3 4 2 5 6];
    case 'H'
        FunctIdx = Map.hFunct;
        ExtractIdx = [Map.hCDR3(2) Map.hGeneName(1) Map.hGeneName(end)];
        SortOrder = [1 2 3];
    case 'L'
        FunctIdx = Map.lFunct;
        ExtractIdx = [Map.lCDR3(2) Map.lGeneName(1) Map.lGeneName(end)];
        SortOrder = [1 2 3];
end

%Delete non-functional genes and keep functional ones
BadLoc = ~strcmpi(VDJdata(:, FunctIdx(1)), 'Y');
for w = 2:length(FunctIdx)
    BadLoc = BadLoc | ~strcmpi(VDJdata(:, FunctIdx(w)), 'Y');
end
BadVDJdata = VDJdata(BadLoc, :);
VDJdata(BadLoc, :) = [];

%Get the cluster cell, which is used to define unique clusters
ClustCell = VDJdata(:, ExtractIdx);
[ClustCell, SortIdx] = sortrows(ClustCell, SortOrder);
VDJdata = VDJdata(SortIdx, :); %Sort by CDR3Lengths mainly
VDJdata(:, Map.GrpNum) = num2cell(zeros(size(VDJdata, 1), 1)); %Reset Groups to 0

%Reduce gene name to family # only
for c = 2:3:size(ClustCell, 2)
    ClustCell(:, c)   = parseGeneName(ClustCell(:, c));   %IGxV#
    ClustCell(:, c+1) = parseGeneName(ClustCell(:, c+1)); %IGxJ#
end

%Convert CDR3 length to string to find unique clusters (matlab workaround)
ClustStr = cell(size(ClustCell, 1), 1); 
for j = 1:size(ClustStr, 1)
    ClustStr{j} = sprintf('%d-%s-%s-', ClustCell{j, :});
end
[~, ~, CrudeClustIdx] = unique(ClustStr, 'stable');

%Determine weighted progress data
CalcPerClust = zeros(1, max(CrudeClustIdx));
for j = 1:max(CrudeClustIdx)
    CalcPerClust(j) = sum(CrudeClustIdx == j)^2;
end
CalcPerClust = cumsum(CalcPerClust / sum(CalcPerClust)) * 100;

%Perform finer clustering per unique crude cluster
GrpNumCt = 0; %buildTreeLink will provide a group number clust from 1 to N
S = buildTreeLink('getheadervar', Map); %For speed, get custom header naming struct S
DelLoc = zeros(size(VDJdata, 1), 1, 'logical'); %Delete duplicate sequences
for j = 1:max(CrudeClustIdx)
    %Perform fine clustering on crude clusters
    ClustLoc = CrudeClustIdx == j;
    Tdata = VDJdata(ClustLoc, :);
    if ~mod(j, 20)
        showStatus(sprintf('  Clustering %d / %d (%0.1f%%)', j, max(CrudeClustIdx), CalcPerClust(j)));
    end

    %If only 1-member, skip everything and update group number only
    if size(Tdata, 1) == 1
        VDJdata{ClustLoc, Map.ChildCount} = 0;
        VDJdata{ClustLoc, Map.GrpNum} = GrpNumCt + 1;
        GrpNumCt = GrpNumCt + 1;
        continue
    end
    
    %Pad sequences CDR3 length also have same Seq Length (required for cluster)
    [Tdata, BadVDJdataT] = padtrimSeqGroup(Tdata, Map, 'cdr3length', 'max', 'Seq');
    if size(BadVDJdataT, 1) > 0
        BadVDJdata = cat(1, BadVDJdata, BadVDJdataT);
    end
    if size(Tdata, 1) == 0 %No valid CDR3, so just mark for deletion and skip
        DelLoc(ClustLoc) = 1;
        continue
    end
    
    %Remove duplicate sequences that can cause issues in clustering
%     Tdata = removeDupSeq(Tdata, VDJheader);
    if size(Tdata, 1) ~= sum(ClustLoc)
        ClustNumIdx = find(ClustLoc);
        DelNumIdx = ClustNumIdx(1:sum(ClustLoc) - size(Tdata, 1)); %Pick how many to delete
        DelLoc(DelNumIdx) = 1;   %mark in VDJdata what to delete
        ClustLoc(DelNumIdx) = 0; %mark in ClustLoc those to keep
    end

    %Find clusters and adjust the group numbers
    [~, Tdata] = buildTreeLink(Tdata, Map, DevPerc, S);  %This will have group numbers 1 to N
    GrpNum2 = cell2mat(Tdata(:, Map.GrpNum)) + GrpNumCt; %This will shift group numbers to unique N+x to N+y
    GrpNumCt = max(GrpNum2);                  %Keep track of highest group num
    Tdata(:, Map.GrpNum) = num2cell(GrpNum2); %Fix the group number
    VDJdata(ClustLoc, :) = Tdata;
end

VDJdata(DelLoc, :) = []; %Delete those that has been reduced from duplicate seq