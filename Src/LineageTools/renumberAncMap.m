%renumberAncMap will renumber AncMap (the ancestry mapping matrix) such
%that the ChildSeqNum (1st col) is sequential and the ParentSeqNum (2nd
%col) is referring back to the correct ChildSeqNum. This is used mainly for
%calculating X-Y coordinates in lineage trees, which can't handle
%out-of-order numbering in column 1.
%
%  AncMap = renumberAncMap(AncMap)
%
%  INPUT
%    AncMap: Mx3 ancestry map of [ChildNum ParentNum Par2ChildDist].
%      Par2ChildDist is either hamming distance or SHM distance from parent
%      to child.
%
%  OUTPUT
%    AncMap: renumbered AncMap based on relative position of child seq num
%      in AncMap.
%
%  NOTE 
%    AncMap does NOT do any resorting, in that the 3rd column remains
%    unchanged. Just purely renumbering 1st and 2nd column.
%
%  EXAMPLE
%    AncMap1 = [
%          1     0     2;
%          3     1     1;
%          5     3     4;
%          7     3     5];
%
%    AncMap2 = renumberAncMap(AncMap1)
%    AncMap2 =
%          1     0     2
%          2     1     1
%          3     2     4
%          4     2     5
%
function AncMap = renumberAncMap(AncMap)
%Adjust parent column (#2) first before child column (#1)
for j = 1:size(AncMap, 1)
    ParIdx = find(AncMap(j, 2) == AncMap(:, 1), 1);
    if ~isempty(ParIdx)
        AncMap(j, 2) = ParIdx;
    else
        AncMap(j, 2) = 0;
    end
end
AncMap(:, 1) = 1:size(AncMap, 1); 