%findTreeClust will return the location of sequence clusters that are
%related by ancestry.
%
%  INPUT
%    AncMap: ancestral map matrix (calcAncMap.m or calcRootedAncMap.m)
%
%  OUTPUT
%    ClustNum: Mx1 matrix where of cluster numbers
%
%  EXAMPLE
%    AncMap = [
%          1     0     1;
%          2     1     2;
%          3     0     1;
%          4     3     1];
%
%    ClustNum = findTreeClust(AncMap)
%    ClustNum =
%          1
%          1
%          2
%          2
%
function ClustNum = findTreeClust(AncMap)
CurNum = 1;
ClustNum = zeros(size(AncMap, 1), 1);
Active = zeros(size(AncMap, 1), 1, 'logical');
Active(1) = 1;
while 1
    ActiveNodes = find(Active); %Initial active node

    %If there are no active nodes, search for the next cluster
    if isempty(ActiveNodes) 
        CurNum = CurNum + 1;
        ActiveNodes = find(ClustNum == 0, 1);
        if isempty(ActiveNodes) %No more to add
            break
        end
    end
    
    %For every active node, search for parent and child
    for q = 1:length(ActiveNodes)
        ClustNum(ActiveNodes(q)) = CurNum;
        ChildLoc  = AncMap(:, 2) == AncMap(ActiveNodes(q), 1);
        ParentLoc = AncMap(:, 1) == AncMap(ActiveNodes(q), 2);
        Active(ChildLoc | ParentLoc) = 1;   %Markup the next nodes
    end
    Active(ClustNum > 0) = 0;  %Prevent rechecking finished checks
end