%calcTreeCoord is used to calculate the node and leaf coordinates
%recursively, given AncMap (the ancestry mapping matrix).
%  
%  TreeCoord = calcTreeCoord(AncMap)
%
%  TreeCoord = calcTreeCoord(AncMap, TreeCoord, ParentNum)  (For recursion)
%
%  INPUT
%    AncMap: Mx3 ancestry map of [ChildNum ParentNum Par2ChildDist].
%      Par2ChildDist is either hamming distance or SHM distance from parent
%      to child. A ParentNum = 0 means germline root.
%   
%  OUTPUT
%    TreeCoord: (M+x)x3 matrix of [Xpos Ypos LeafStatus].
%      IF the distance from root to 1st real sequence is 0, then x = 0.
%      IF the root is not the same as th 1st sequence, then x = 1 and a
%        extra coordinate is provided from root to 1st sequence.
%      LeafStatus is 0 if this child is also a parent.
%      LeafStatus is 1 if this child ends the tree as a leaf.
%
%  NOTE
%    This uses recursion to fill in TreeCoord, hence TreeCoord is a hidden
%    input to the function that shouldn't be used by the user.

function TreeCoord = calcTreeCoord(AncMap, varargin)
if nargin == 1 %Initial summon
    if any(findTreeCycle(AncMap))
        error('%s: Cannot have cyclic dependency in AncMap.', mfilename);
    end
    AncMap = renumberAncMap(AncMap);
    TreeCoord = zeros(size(AncMap, 1), 3);
    ParentNum = max([AncMap(1, 2), 1]);
elseif nargin == 3 %Recursive summon
    TreeCoord = varargin{1};
    ParentNum = varargin{2};
else
    error('%s: Incorrect number of inputs.', mfilename);
end

%If it is a leaf (no child), update TreeCoord. If not, resummon until you do.
ChildNum = findChild(AncMap, ParentNum);
if isempty(ChildNum) %Leaf
    TreeCoord(ParentNum, 2) = max(TreeCoord(:, 2)) + 1; %Count the Y value
    TreeCoord(ParentNum, 3) = 1; %Mark as leaf;
else %Node
    CurMaxY = max(TreeCoord(:, 2)) + 1; %The Y value that you WOULD have accepted, but must now wait till all childs are added.
    for j = 1:length(ChildNum)
        TreeCoord(ChildNum(j), 1) = TreeCoord(ParentNum, 1) + AncMap(ChildNum(j), 3); %Update the X value relative to parent's position
        TreeCoord = calcTreeCoord(AncMap, TreeCoord, ChildNum(j));                    %Find the child's child to update TreeCoord
    end
    NewMaxY = max(TreeCoord(:, 2));                     %This is the new maximum Y position, after all childen were processed.
    TreeCoord(ParentNum, 2) = (CurMaxY + NewMaxY) / 2;  %Now take the mid point after all children were processed.
end