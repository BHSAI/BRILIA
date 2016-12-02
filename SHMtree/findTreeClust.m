function ClustMap = findTreeClust(AncMap,Option)
% if Option == 1 %Find cluster by connection

%Build the starting cluster
ClustNum = 1;
ClustMap = [[1:size(AncMap,1)]' zeros(size(AncMap,1),2)]; %[SeqNumber ClusterNum CurrentNodeAssessed]
%ClustMap(AncMap(:,2)==0,2) = ClustNum;
ClustMap(1,3) = 1; %Start with an active node

k = 1; %Kill counter to prevent inf loops
while k <= size(ClustMap,1)
    ActiveNode = find(ClustMap(:,3) == 1); %Initial active node
    
    %If there are no active nodes, search for the next one
    if isempty(ActiveNode) 
        ClustNum = ClustNum + 1;
        NextLoc = find(ClustMap(:,2)==0); %Search for unasigned seq
        if isempty(NextLoc) %No more nodes, so just end
            break
        end
        ClustMap(NextLoc(1),3) = 1; %Mark the next node
        continue
    end
    
    %For every active node, search for parent and child
    for q = 1:length(ActiveNode)
        ClustMap(ActiveNode(q),2) = ClustNum; %Save cluster num
        ClustMap(ActiveNode(q),3) = 0; %Deactivate current node

        SeqIsParentOf = find(AncMap(:,2) == AncMap(ActiveNode(q),1)); %Find child of current seq
        SeqIsChildOf = find(AncMap(:,1) == AncMap(ActiveNode(q),2)); %Find parent of current seq
        
        ClustMap([SeqIsParentOf; SeqIsChildOf],3) = 1; %Markup the next node        
    end
    ClustMap((ClustMap(:,2)>0),3) = 0; %Prevents working on already-worked nodes.
    
    k = k+1;
end
ClustMap = ClustMap(:,1:2);