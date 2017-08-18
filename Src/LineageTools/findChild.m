%findChild will return the numbers for the children of sequence.

function ChildLoc = findChild(AncMap,ParentNum)
ChildIdx = AncMap(:,2) == ParentNum;
ChildLoc = AncMap(ChildIdx,1);
