%findTreeCycle will look for the AncMap entries that causes a cyclic loop.
%
%  CycleLoc = findTreeCycle(AncMap)

function CycleLoc = findTreeCycle(AncMap)
%To find a cycle, essentially remove leaves until no more is left. What's
%left is the cyclic dependencies.


CycleLoc = ones(size(AncMap,1),1) == 1;
SumLeaf = sum(CycleLoc);
while 1
    for j = 1:size(AncMap,1)
        if CycleLoc(j) == 1
            ChildLoc = AncMap(:,2) == AncMap(j,1);
            if max(ChildLoc) == 0
                AncMap(j,2) = -AncMap(j,2);
                CycleLoc(j) = 0;
            end
        end
    end
    if sum(CycleLoc) == SumLeaf
        break
    else
        SumLeaf = sum(CycleLoc);
    end
end
% 
%     CurLeaf = sum(CycleLoc);
% 
% 
% 
% 
% 
% 
% %Remove those that lack a child
% 
% %Evaluate those are internal
% 
% 
% 
% Laps = zeros(size(AncMap,1),1);
% for k = 1:2
%     for j = 1:size(AncMap,1)
%         ParentLoc = AncMap(:,1) == AncMap(j,2);
%         if max(ChildLoc) > 0 && max(ParentLoc) > 0
%             Laps(j) = Laps(j)+1;
%         end
%     end
% end
% 
% CycleLoc = Laps >= 2;
% % 
% % %Check for cyclic dependencies, and end the cycle now
% % for j = 1:length(ChildLocT)
% %     if CycleLoc(ChildLocT(j)) == 1
% %         AncMap(ChildLocT(j),2) = 0;
% %         CyclicCut = 1
% %     end
% % end
% % CycleLoc(ChildLocT) = 1;
% % 
% % for j = 1:length(ChildLocT)
% %     [CycleLoc, CyclicCut] = findTreeChild(AncMap,ChildLocT(j),CyclicCut,CycleLoc);
% % end
% % 
% % if nargout == 2
% %     varargout{1} = CyclicCut;
% % end
