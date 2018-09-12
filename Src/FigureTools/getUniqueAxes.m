%getUniqueAxes is designed to find unique axes, which can be cumbersome if
%using a workaround where one plot has > 1 axes overlaying its position. If
%two axes have the same positions, getUniqueAxes will return the one unique
%axes, but for each one return which are supposedly the duplicate one.
%Returns the position and tight insets of each unique plot as if the axes
%formed one plot.
%
%  INPUT
%    Gx: figure handle
% 
%  OUTPUT
%    UnqAxCell: a cell of axes handle array, where each cell stores axes
%      handles that have > 50% overlap in plotting area.
%    AllPositions: Mx4 position matrix of the first plot with a background
%      color per each unique group of axes
%    AllTights: Mx4 tight inset matrix of the first plot with a background
%      color per each unique group of axes

function [UnqAxCell, AllPositions, AllTights] = getUniqueAxes(Gx)
if ~(ishandle(Gx) && strcmpi(class(Gx), 'matlab.ui.Figure'))
   error('%s: Must input a figure handle', mfilename);
end
Axs = findobj(Gx, 'type', 'axes');

%Determine locations of all axes and tight insets
AllTights  = zeros(length(Axs), 4);
AllPositions = zeros(length(Axs), 4);
AllColor = ones(length(Axs), 1); %1 if color is NOT 'none', and 0 if it is
for a = 1:length(Axs)
    AllTights(a, :) = get(Axs(a), 'TightInset');
    AllPositions(a, :) = get(Axs(a), 'Position');
    AxColor = get(Axs(a), 'Color');
    if ischar(AxColor) && strcmpi(AxColor, 'none')
        AllColor(a) = 0;
    end
end

%Determine all unique positions when using a workaround overlay axes.
UnqAxNum = zeros(length(Axs), 1);
u = 1; %Unique counter
for j = 1:length(Axs)
    if UnqAxNum(j) > 0; continue; end %Already assigned a unique num
    UnqAxNum(j) = u;
    Box1 = AllPositions(j, :);
    Box1Area = Box1(3) * Box1(4);
    for k = 2:length(Axs)
        if UnqAxNum(k) > 0; continue; end %Already assigned a unique num
        Box2 = AllPositions(k, :);
        Box2Area = Box2(3) * Box2(4);
        Ho = min(Box1(2) + Box1(4), Box2(2) + Box2(4)) - max(Box1(2), Box2(2));
        Wo = min(Box1(1) + Box1(3), Box2(1) + Box2(3)) - max(Box1(1), Box2(1));
        if Ho < 0 && Wo < 0; continue; end %nothing overlaps
        if (Ho*Wo) / max(Box2Area, Box1Area) > 0.50 %Consider this the same
            UnqAxNum(k) = u;
        end
    end
    u = u + 1;
end

%Get all the unique Ax , Positions, and Tight Insets
KeepIdx = zeros(length(Axs), 1, 'logical');
UnqAxCell = cell(max(UnqAxNum), 1);
for j = 1:max(UnqAxNum)
    Idx = find(UnqAxNum == j);
    KeepIdx(Idx(1)) = 1;
    UnqAxCell{j} = Axs(Idx);
    AllTights(Idx(1), :) = max(AllTights(Idx, :), [], 1);
end
AllPositions = AllPositions(KeepIdx, :);
AllTights = AllTights(KeepIdx, :);
