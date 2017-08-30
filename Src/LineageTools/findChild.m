%findChild will return the locations numbers for the children of sequence.

function ChildLoc = findChild(AncMap, ParentNum)
ChildLoc = AncMap(AncMap(:,2) == ParentNum, 1);
