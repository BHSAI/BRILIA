%clusterByJunction will cluster VDJ sequences by the V+J gene family or
%number, and the CDR3 lengths. This only assigns a cluster by number, but
%lineage-based clustering must be done immediately afterward.
%
%  [VDJdata, BadVDJdata] = clusterByJunction(VDJdata, Map, varargin)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    Map: header map of BRILIA data cells
%    'fam': family number. Empty is default for up to gene number.
%    
%  OUTPUT
%    VDJdata: main BRILIA data cell with the GrpNum filled by V+J+CDR3Len
%
%  NOTE
%    Default clustering is by Gene Number such as 'IGHV1-11' and 'IGHJ1'.
%    Use 'fam' to cluster by 'IGHV1' and 'IGJ1'.

function [VDJdata, BadVDJdata] = clusterByJunction(VDJdata, Map, varargin)
UseFamNum = startsWith('fam', varargin, 'ignorecase', true);

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
        if UseFamNum %Use the family number, ex: IGHV1
            ClustCell(:, c) = parseGeneName(ClustCell(:, c));
        else %Default. Use the gene number, ex: IGHV1-11
            [FamNum, GeneNum] = parseGeneName(ClustCell(:, c));
            ClustCell(:, c) = strcat(FamNum, GeneNum);
        end
        UnqStrPat{c} = '%s-';
    end
end
UnqStrPat = [UnqStrPat{:}];
CrudeClustStr = cell(size(ClustCell, 1), 1); %matlab workaround for unique rows of cells
for j = 1:size(CrudeClustStr, 1)
    CrudeClustStr{j} = sprintf(UnqStrPat, ClustCell{j, :});
end
[~, ~, CrudeClustIdx] = unique(CrudeClustStr, 'stable');
[CrudeClustIdx, SortIdx] = sort(CrudeClustIdx);
VDJdata = VDJdata(SortIdx, :);
VDJdata(:, Map.GrpNum) = num2cell(CrudeClustIdx);
