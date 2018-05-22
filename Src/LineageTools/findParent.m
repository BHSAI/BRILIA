%findParent will return the location numbers for the parent of sequence.
%
%  ParentNum = findParent(AncMap, ChildNum)
%
%  INPUT
%    AncMap: ancestry map
%    ChildNum: child number for which the parent number is sought
%
%  OUTPUT
%    ParentNum: the parent number(s) of the child

function ParentNum = findParent(AncMap, ChildNum)
ParentNum = AncMap(AncMap(:,1) == ChildNum, 2);
