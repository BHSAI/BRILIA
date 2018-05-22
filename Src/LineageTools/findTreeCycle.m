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
SumLeaf = length(CycleLoc);
while 1
    %Remove those that are leaves
    for j = 1:size(AncMap, 1)
        if CycleLoc(j)
            ChildLoc = AncMap(:, 2) == AncMap(j, 1);
            if max(ChildLoc) == 0
                AncMap(j, 2) = -AncMap(j, 2); %Set negative to prevent redoing search here
                CycleLoc(j) = 0;
            end
        end
    end
    
    %Keep track of current and previous leaf count to determine loop break
    if sum(CycleLoc) == SumLeaf
        break
    end
    SumLeaf = sum(CycleLoc);
end