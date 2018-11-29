%findChild will return the locations numbers for the children of sequence.
%
%  ChildNum = findChild(AncMap, ParentNum)
%
%  INPUT
%    AncMap: ancestry map
%    ParentNum: parent number for which the child(ren) number(s) is(are) 
%      sought 
%
%  OUTPUT
%    ChildNum: the child number(s) of the parent

function ChildNum = findChild(AncMap, ParentNum)
ChildNum = AncMap(AncMap(:,2) == ParentNum, 1);
