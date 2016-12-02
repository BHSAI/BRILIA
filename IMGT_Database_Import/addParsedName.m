
%Properly sort the family names
function Xmap = addParsedName(Xmap)

AddParsedName = cell(size(Xmap,1),4);
for j = 1:size(Xmap,1)
    AddParsedName(j,:) = deal(parseGeneName(Xmap{j,3}));
end
Xmap = cat(2,Xmap,AddParsedName);
% 
% Xmap = sortrows(Xmap,4);
end
