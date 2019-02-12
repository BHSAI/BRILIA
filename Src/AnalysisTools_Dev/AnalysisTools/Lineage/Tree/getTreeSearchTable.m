%getTreeSearchTable is designed to display a table listing the tree size,
%V, D, J family, and tree image file name for use as a quick glance of the
%lineage trees for a particular file.
%
%  TreeSearchTable = getTreeSearchTable(VDJdata, VDJheader)
%
%  INPUT
%    TreeData: output from collectTreeData
%
%  OUTPUT
%    TreeSearchTable: condensed table for viewing contents of the tree. The
%      first row is the header of the column.
%        GrpNum: BRILIA-assigned group number
%        TreeSize: total number of unique sequences
%        TotalTemplateCount: sum of all sequence template counts
%        H-Vgene: heavy chain V gene name
%        H-Dgene: heavy chain D gene name
%        H-Jgene: heavy chain J gene name
%        H-CDR3: germline CDR3 per cluster of the heavy chain
%        L-Vgene: light chain V gene name
%        L-Jgene: light chain J gene name
%        L-CDR3: germline CDR3 per cluster of the light chain
%        TreeImageName: file name of the corresponding group lineage tree (to be filled later)

function TreeSearchTable = getTreeSearchTable(VDJdata, VDJheader)
%Determine what data to extract from VDJdata
Map = getVDJmapper(VDJheader);
if strcmpi(Map.Chain, 'H')
    ExtractThese = [Map.GrpNum; Map.Template; Map.Template; Map.hCDR3(1);  Map.hGeneName]; %The 2 Map.Template is a workaround. 2nd column will be replaced.
elseif strcmpi(Map.Chain, 'L')
    ExtractThese = [Map.GrpNum; Map.Template; Map.Template; Map.lCDR3(1);  Map.lGeneName];
else
    ExtractThese = [Map.GrpNum; Map.Template; Map.Template; Map.hCDR3(1);  Map.hGeneName(1); Map.lCDR3(1); Map.lGeneName];
end

%For each cluster, collect data
GrpNum = cell2mat(VDJdata(:, Map.GrpNum));
UnqGrpNum = unique(GrpNum);
TreeSearchTable = cell(length(UnqGrpNum), length(ExtractThese) + 1);
for y = 1:length(UnqGrpNum)
    GrpIdx = find(UnqGrpNum(y) == GrpNum);
    TreeSearchTable(y, 1:end-1) = VDJdata(GrpIdx(1), ExtractThese);
    TreeSearchTable{y, 2} = length(GrpIdx);
    TreeSearchTable{y, 3} = sum(cell2mat(VDJdata(GrpIdx, Map.Template)));
end
TreeSearchHeader = [VDJheader(ExtractThese) 'TreeImageName'];
TreeSearchHeader{2} = 'UniqueSequences';

TreeSearchTable = cat(1, TreeSearchHeader, TreeSearchTable);
