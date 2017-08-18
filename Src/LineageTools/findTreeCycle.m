%findTreeCycle will look for the AncMap entries that causes a cyclic loop.
%
%  CycleLoc = findTreeCycle(AncMap)

function CycleLoc = findTreeCycle(AncMap)
%To find a cycle, essentially remove leaves until no more is left. What's
%left is the cyclic dependencies.
CycleLoc = ones(size(AncMap,1),1) == 1;
SumLeaf = sum(CycleLoc);
while 1
    for j = 1:size(AncMap,1)
        if CycleLoc(j) == 1
            ChildLoc = AncMap(:,2) == AncMap(j,1);
            if max(ChildLoc) == 0
                AncMap(j,2) = -AncMap(j,2);
                CycleLoc(j) = 0;
            end
        end
    end
    if sum(CycleLoc) == SumLeaf
        break
    else
        SumLeaf = sum(CycleLoc);
    end
end
