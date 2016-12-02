%calcTreeCoord is used to calculate the leaf coordinates, using a while
%within a while loop method. 

function TreeCoord = calcTreeCoord(AncMap,varargin)
if isempty(varargin)
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

%1) Find the children
ChildNum = find(ParentNum == AncMap(:,2));

%2) If it is a leaf, update TreeCoord. If not, resummon until you do.
if isempty(ChildNum)
    TreeCoord(ParentNum,2) = max(TreeCoord(:,2)) + 1; %Count the Y value
    TreeCoord(ParentNum,3) = 1; %Mark as leaf;
else
    CurMaxY = max(TreeCoord(:,2))+1;
    for j = 1:length(ChildNum)
        TreeCoord(ChildNum(j),1) = TreeCoord(ParentNum,1) + AncMap(ChildNum(j),3); %Count the X value     
        TreeCoord = calcTreeCoord(AncMap,TreeCoord,ChildNum(j));
    end
    NewMaxY = max(TreeCoord(:,2));
    TreeCoord(ParentNum,2) = (CurMaxY + NewMaxY) / 2; %Count the Y value
end