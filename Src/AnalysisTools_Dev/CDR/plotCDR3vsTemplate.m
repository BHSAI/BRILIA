%plotCDR3vsTemplate will plot the number of unique CDR3 per clonotype vs
%the total template count within that clonotype.

function [M, UnqLen] = plotCDR3vsTemplate(Data)
Len = arrayfun(@(x) length(x.HCDR3{1}), Data)';
KeepLoc = Len > 1;

X = [Data.Template]';
Y = [Data.HCDR3Count]';
Len = arrayfun(@(x) length(x.HCDR3{1}), Data)';
UnqLen = unique(Len);

ClrMap = jet*0.75;
ClrIdx = 1:floor((size(ClrMap, 1)-1)/length(UnqLen)):size(ClrMap, 1);

M = zeros(length(UnqLen), 1);
for j = 1:length(UnqLen)
    Idx = Len == UnqLen(j) & KeepLoc;
    scatter(X(Idx), Y(Idx), 20, ClrMap(ClrIdx(j), :), 'fill');
    hold on
    M(j) = X(Idx)\Y(Idx);
end
hold off

