%findChild will return the numbers for the parent of sequence.

function ParentLoc = findParent(AncMap,ChildNum)
ParentIdx = AncMap(:,1) == ChildNum;
ParentLoc = AncMap(ParentIdx,2);
