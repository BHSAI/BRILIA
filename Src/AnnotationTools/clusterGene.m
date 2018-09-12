%clusterGene will do the full clustering of sequences based on lineage
%viability, distances, germline sequence annotation, etc.
%
%PROCEDURE
%1) Group crude clonotypes by unique CDR3Length, Vfamily, and Jfamily
%   For both H and L chains.
%2) Compute SHM pairwise sequence distances and then use the Motif and Mut
%   scores to "block" certain lineage linkings. Essentially, having more
%   coldspot than hotspot mutations raises a red flag for linking.
%3) Link by nearest distance first to get the fine clusters.
%4) Within each cluster, look for the closest SHM distance to the germline.
%     Closest Distance >> 
%     Higher Template Count >> 
%     Smaller Distance to others members >>
%5) Determine if clusters can be linked to the same germline 
%
%
%  VDJdata = clusterGene(VDJdata, Map)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    Map: header map of BRILIA data cells
%
%  OUTPUT
%    VDJdata: modified VDJdata where sequences are clustered based on
%      lineage relationships.

function [VDJdata, BadVDJdata] = clusterGene(VDJdata, Map)
error('%s: this is obsolete. See clusterGeneByLineage', mfilename);

%Determine chain and extract key locations, sorting in order
switch Map.Chain
    case 'HL'
        FunctIdx = [Map.hFunct; Map.lFunct];
        CDR3LIdx = [Map.hCDR3(2); Map.lCDR3(2)];
        GeneNIdx = [Map.hGeneName([1 end]); Map.lGeneName([1 end])];
    case 'H'
        FunctIdx = Map.hFunct;
        CDR3LIdx = Map.hCDR3(2);
        GeneNIdx = Map.hGeneName([1 end]);
    case 'L'
        FunctIdx = Map.lFunct;
        CDR3LIdx = Map.lCDR3(2);
        GeneNIdx = Map.lGeneName([1 end]);
end

%Delete non-functional genes and keep functional ones
BadLoc = any(~strcmpi(VDJdata(:, FunctIdx), 'Y'), 2) | any(cellfun('isempty', VDJdata(:, CDR3LIdx)), 2);
BadVDJdata = VDJdata(BadLoc, :);
VDJdata = VDJdata(~BadLoc, :);

%Get the cluster cell, which is used to define unique clusters
ClustCell = VDJdata(:, [CDR3LIdx; GeneNIdx]);
[ClustCell, SortIdx] = sortrows(ClustCell);
VDJdata = VDJdata(SortIdx, :); %Sort by CDR3Lengths mainly
VDJdata(:, Map.GrpNum) = {0};
VDJdata(:, Map.ChildCount) = {0};
VDJdata(:, Map.ParNum) = {0};

%Reduce gene name to family # only, and get unique clusters
UnqStrPat = repelem({'%03d-'}, 1, size(ClustCell, 2));
for c = 1:size(ClustCell, 2)
    if ~isnumeric(ClustCell{1, c})
        ClustCell(:, c) = parseGeneName(ClustCell(:, c));
        UnqStrPat{c} = '%s-';
    end
end
UnqStrPat = [UnqStrPat{:}];
CrudeClustStr = cell(size(ClustCell, 1), 1); %matlab workaround for unique rows of cells
for j = 1:size(CrudeClustStr, 1)
    CrudeClustStr{j} = sprintf(UnqStrPat, ClustCell{j, :});
end
[~, ~, CrudeClustIdx] = unique(CrudeClustStr, 'stable');

%Calculate code progress % increments
Progress = zeros(1, max(CrudeClustIdx));
for j = 1:max(CrudeClustIdx)
    Progress(j) = sum(CrudeClustIdx == j)^2; %Scales roughly with N^2
end
Progress = cumsum(Progress) / sum(Progress) * 100;

%Perform finer clustering per unique crude cluster
GrpNumStart = 1; %buildTreeLink will provide a group number clust from 1 to N
DelLoc = zeros(size(VDJdata, 1), 1, 'logical'); %Delete sequence with unreasonable CDR3s (this should be pushed elsewhere...)
for j = 1:max(CrudeClustIdx)
    if ~mod(j, 20)
        showStatus(sprintf('  Clustering %d / %d (%0.1f %%)', j, max(CrudeClustIdx), Progress(j)));
    end
    
    %If only 1-member, skip everything and update group number only
    ClustLoc = CrudeClustIdx == j;
    if sum(ClustLoc) == 1
        VDJdata{ClustLoc, Map.GrpNum} = GrpNumStart;
        GrpNumStart = GrpNumStart + 1;
        continue
    end
    
    %Pad sequences CDR3 length also have same Seq Length (required for cluster)
    Tdata = VDJdata(ClustLoc, :);
    [Tdata, BadVDJdataT] = padtrimSeqGroup(Tdata, Map, 'cdr3length', 'max', 'Seq');
    if ~isempty(BadVDJdataT)
        BadVDJdata = cat(1, BadVDJdata, BadVDJdataT);
    end
    if isempty(Tdata) %No valid CDR3. Mark for deletion and skip.
        DelLoc(ClustLoc) = 1;
        continue
    end
    
    if size(Tdata, 1) ~= sum(ClustLoc)
        ClustNumIdx = find(ClustLoc);
        DelCount = sum(ClustLoc) - size(Tdata, 1);
        DelLoc(ClustNumIdx(1:DelCount)) = 1;   %mark in VDJdata what to delete
        ClustLoc(ClustNumIdx(1:DelCount)) = 0; %adjust the "save to" locations
    end

    [VDJdata(ClustLoc, :), GrpNumStart] = buildTreeLink(Tdata, Map, GrpNumStart);
end

VDJdata(DelLoc, :) = []; %Delete those that has been reduced from duplicate seq

%Sort by group number and parent number
[~, SortIdx] = sortrows(cell2mat(VDJdata(:, [Map.GrpNum Map.ParNum])));
VDJdata = VDJdata(SortIdx, :);
