%findTreeCycle will look for the AncMap entries that causes a cyclic loop.
%
%  CycleLoc = findTreeCycle(AncMap)
%
%  INPUT
%    AncMap: ancestral map index from calcAncMap
%
%  OUTPUT
%    CycleLoc: logical array of size(AncMap, 1)x1 marking members that have
%      parent-child-parent cyclic dependencies
%
%  EXAMPLE
%    AncMap = [
%          1     4     1;
%          2     1     1;
%          3     2     1;
%          4     3     1];
%
%    CycleLoc = findTreeCycle(AncMap)
%    CycleLoc =
%          1
%          1
%          1
%          1

function CycleLoc = findTreeCycle(AncMap)
%To find a cycle, remove leaves until no more is left.
CycleLoc = ones(size(AncMap, 1), 1, 'logical');
CurSumLeaf = length(CycleLoc);
SumLeaf = 0;
while SumLeaf ~= CurSumLeaf
    %Remove those that are leaves
    for j = 1:size(AncMap, 1)
        if CycleLoc(j)
            ChildLoc = AncMap(:, 2) == AncMap(j, 1);
            if ~any(ChildLoc)
                AncMap(j, 2) = -AncMap(j, 2); %Set negative to prevent redoing search here
                CycleLoc(j) = 0;
            end
        end
    end
    
    SumLeaf = CurSumLeaf;
    CurSumLeaf = sum(CycleLoc);
end