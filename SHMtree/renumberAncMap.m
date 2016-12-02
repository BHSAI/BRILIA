%renumberAncMap will get an AncMap (ancestram mapping matrix), and then
%renumber them so that the ChildSeq (1st col) goes from 0 to 1, and the
%child's parent loc matches, relatively.

%  EX: AncMap =  [1 0;  3 1;  5 3;  7 3];
% 
%      AncMap = renumberAncMap(AncMap)
%        AncMap =
%              1     0
%              2     1
%              3     2
%              4     2

function AncMap = renumberAncMap(AncMap)

%Make sure that AncMap is numbered 1 to N on first column.
for j = 1:size(AncMap,1)
    ParRelLoc = find(AncMap(j,2) == AncMap(:,1));
    if ~isempty(ParRelLoc )
        AncMap(j,2) = ParRelLoc;
    else
        AncMap(j,2) = 0;
    end
end
AncMap(:,1) = 1:size(AncMap,1);
