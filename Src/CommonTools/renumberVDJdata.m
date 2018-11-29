%renumberVDJdata will renumber AND reorder the groups and/or sequences. 
%
%  VDJdata = renumberVDJdata(VDJdata, Map, Option, Option)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    Map: structure of BRILIA data header indices
%    Option ['grp']: Renumber group numbers. (Default if no Option given)
%    Option ['seq']: Renumber sequence numbers. 
%
%  OUTPUT
%    VDJdata: modified VDJdata with renumbered groups and/or sequences. 
%      Data will be sorted by group number, but will preserve sequence
%      number order since legacy BRILIA versions relies on 1st sequence of
%      a group being the root of a tree.
%
%  EXAMPLES
%    VDJdata = num2cell([[1 3 4 7 2 8 13]' [1 1 1 3 1 3 3]']);
%    VDJheader = {'SeqNum', 'GroupNum'};
%    Map = getVDJmapper(VDJheader);
%
%    VDJdata = renumberVDJdata(VDJdata, Map, 'seq')
%    SeqNum    GrpNum  -->   SeqNum    GrpNum
%    1         1             1         1
%    3         1             2         1
%    4         1             3         1
%    7         3             4         1
%    2         1             5         3
%    8         3             6         3
%    13        3             7         3
%
%    VDJdata = renumberVDJdata(VDJdata, Map, 'grp')
%    SeqNum    GrpNum  -->   SeqNum    GrpNum
%    1         1             1         1
%    3         1             3         1
%    4         1             4         1
%    7         3             2         1
%    2         1             7         2
%    8         3             8         2
%    13        3             13        2
%
%    VDJdata = renumberVDJdata(VDJdata, Map, 'seq', 'grp')
%    SeqNum    GrpNum  -->   SeqNum    GrpNum
%    1         1             1         1
%    3         1             2         1
%    4         1             3         1
%    7         3             4         1
%    2         1             5         2
%    8         3             6         2
%    13        3             7         2

function VDJdata = renumberVDJdata(VDJdata, Map, varargin)
if iscell(Map)
    Map = getVDJmapper(Map);
end
if isempty(varargin)
    RenumSeq = 0;
    RenumGrp = 1;
else
    RenumSeq = contains('seq', varargin, 'IgnoreCase', true);
    RenumGrp = contains('grp', varargin, 'IgnoreCase', true);
end

GrpNum = cell2mat(VDJdata(:, Map.GrpNum));
if RenumGrp
    [~, ~, GrpNum] = unique(GrpNum, 'stable');
    VDJdata(:, Map.GrpNum) = num2cell(GrpNum);
end
[~, SortIdx] = sort(GrpNum);
VDJdata = VDJdata(SortIdx, :);

if RenumSeq
    VDJdata(:, Map.SeqNum) = num2cell(1:size(VDJdata, 1));
end