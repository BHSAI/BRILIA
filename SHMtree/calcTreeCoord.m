%calcTreeCoord is used to calculate the node and leaf coordinates
%recursively, given AncMap (the ancestry mapping matrix).
%  
%  TreeCoord = calcTreeCoord(AncMap)
%
%  TreeCoord = calcTreeCoord(AncMap,TreeCoord,ParentNum)  (For recursion)
%
%  INPUT
%    AncMap: Mx3 ancestry map of [ChildNum ParentNum Par2ChildDist].
%      Par2ChildDist is either hamming distance or SHM distance from parent
%      to child.
%   
%  OUTPUT
%    TreeCoord: Mx3 matrix of [Xpos Ypos NodeStatus]. LeafStatus is 0 if
%      this child is also a parent, and Nodestatus is 1 if this child ends
%      the tree as a leaf.
%
%  NOTE
%    This uses recusion to fill in TreeCoord, hence TreeCoord is a hidden
%      input to the function that shouldn't be used by the user.

function TreeCoord = calcTreeCoord(AncMap,varargin)

%Determine if it's the initial or recursive function summon
if isempty(varargin)
    %Renumber AncMap the first time to ensure ChildNum is sequential
    AncMap = renumberAncMap(AncMap);
    TreeCoord = zeros(size(AncMap,1),3); %Xcoor, Ycoor, LeafStatus (1 or 0)
    ParentNum = AncMap(1,2);
    if ParentNum == 0
        ParentNum = ParentNum+1;
    end
else
    TreeCoord = varargin{1};
    if length(varargin) == 2
        ParentNum = varargin{2};
    end
end

%Find the children of this parent
ChildNum = find(ParentNum == AncMap(:,2));

%If it is a leaf (no child), update TreeCoord. If not, resummon until you do.
if isempty(ChildNum) %Leaf
    TreeCoord(ParentNum,2) = max(TreeCoord(:,2)) + 1; %Count the Y value
    TreeCoord(ParentNum,3) = 1; %Mark as leaf;
else %Node
    CurMaxY = max(TreeCoord(:,2)) + 1; %The Y value that you WOULD have accepted, but must now wait till all childs are added.
    for j = 1:length(ChildNum)
        TreeCoord(ChildNum(j),1) = TreeCoord(ParentNum,1) + AncMap(ChildNum(j),3); %Update the X value relative to parent's position
        TreeCoord = calcTreeCoord(AncMap,TreeCoord,ChildNum(j)); %Find the child's child to update TreeCoord
    end
    NewMaxY = max(TreeCoord(:,2)); %This is the new maximum Y position, after all childen were processed.
    TreeCoord(ParentNum,2) = (CurMaxY + NewMaxY) / 2; %Now take the mid point after all children wer processed.
end