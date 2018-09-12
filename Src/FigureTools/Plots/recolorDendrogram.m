%recolorDendrogram will recolor the dendrogram according to a colormap and
%grouping number.
%
%  recolorDendrogram(Dx, Grouping, ColorMap)
%
%  INPUT
%    Dx: handle returned from the dendrogram ( ex: Dx = dendrogram(M) )
%    Grouping: Sequential order of leaves but assigned to a group. ( ex:
%      [1 1 1 3 3 3 2 3 3 3 2] )
%    ColorMap: Mx3 color map matrix

function recolorDendrogram(Dx, Grouping, ColorMap)
if nargin < 3
    ColorMap = jet;
end
UnqGroup = unique(Grouping);

%Color coding by group number
MapIncr = floor((size(ColorMap, 1)-1) / length(UnqGroup));
XLabel = get(gca, 'XTickLabel');
XLabelNum = zeros(size(XLabel, 1), 1);
for j = 1:length(XLabelNum)
    XLabelNum(j) = convStr2Num(XLabel(j, :));
end

%Generate the color code in order of the dendrogram leaves
ClrCode = zeros(length(XLabelNum), 3);
for j = 1:length(UnqGroup)
    ArmPosNum = find(Grouping == UnqGroup(j));
    [~, XPosNum] = intersect(XLabelNum, ArmPosNum);
    ClrCode(XPosNum, :) = repmat(ColorMap(1+MapIncr*(j-1), :), length(XPosNum), 1);
end

for k = 1:length(Dx)
    Dx(k).Color = ClrCode(k, :);
end