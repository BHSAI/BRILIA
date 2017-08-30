%findChild will return the location numbers for the parent of sequence.

function ParentLoc = findParent(AncMap, ChildNum)
ParentLoc = AncMap(AncMap(:,1) == ChildNum, 2);
