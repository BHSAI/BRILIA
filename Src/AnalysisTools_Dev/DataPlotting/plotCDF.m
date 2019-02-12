function plotCDF(O, varargin)
if ~isfield(O.Ax, 'CDF') || isempty(O.Ax.CDF) || endsWith(class(O.Ax.CDF(1)), 'GraphicsPlaceholder') || ~isvalid(O.Ax.CDF(1))
    figure;
    O.Ax.CDF = axes;
end

N = numel(O.ModData);
StdColor = getStdColor;

% LxCDF = gobjects(N, 1);
% LxYIntercept = gobjects(N, 1); 
CDFs = zeros(N, 1);


Cutoff = 5;

for j = 1:N
    Template = sort(O.ModData(j).Data);
    CDF = (1:numel(Template))/numel(Template);
    plot(O.Ax.CDF, Template, CDF, 'Color', StdColor(O.ModData(j).Group, :), 'LineWidth', 1.5);
    if j == 1
        hold(gca, 'on')
    end
    
    %Computer the y-intercept line
    YIdx = find(Template <= Cutoff, 1, 'last');
    CDFs(j) = CDF(YIdx);
%    LxYIntercept(j) = scatter(O.Ax.CDF, Template(YIdx), CDFs(j) , 'o', 'MarkerEdgeColor', StdColor(O.Group(j), :));
end
% plot([Cutoff Cutoff], [0 1], 'k--');
hold(gca, 'off')
xlim([0 50])

Group = [O.ModData.Group]';

Mean = splitapply(@mean, CDFs, Group);
Std = splitapply(@std, CDFs, Group);
Entropy = zeros(N, 1);
for j = 1:numel(Entropy)
    [~, Entropy(j)] = calcDiversity('Shannon', O.ModData(j).Data);
end
AvgEntropy = splitapply(@mean, Entropy, Group);
StdEntropy = splitapply(@std, Entropy, Group);

TotalS = cellfun('length', {O.ModData.Data})';
Evenness = Entropy ./ log(TotalS);
AvgEvenness = splitapply(@mean, Evenness, Group);
StdEvenness = splitapply(@std, Evenness, Group);

%Assemble text legend
DataHeader = {'Group' sprintf('Fr. %s Cutoff (%d)', char(8804), Cutoff) 'Shannon Entropy', 'Evenness'};

%Final Adjustements
setPlotTickDecimal(O.Ax.CDF, 0, 1);
title(O.Ax.CDF, [O.DataName 'CDF']);
resizeSubplots(O.Ax.CDF);
return 

GroupTxt = cellfun(@(x) sprintf('%d', x), num2cell(unique(O.Group)), 'un', 0);
CutoffTxt = cellfun(@(x,y) sprintf('%0.3f %s %0.3f', x, char(177), y), num2cell(Mean), num2cell(Std), 'un', 0);
ShannonTxt = cellfun(@(x,y) sprintf('%0.1f %s %0.1f', x, char(177), y), num2cell(AvgEntropy), num2cell(StdEntropy), 'un', 0);
EvennessTxt = cellfun(@(x,y) sprintf('%0.3f %s %0.3f', x, char(177), y), num2cell(AvgEvenness), num2cell(StdEvenness), 'un', 0);
DataValue = [GroupTxt, CutoffTxt, ShannonTxt, EvennessTxt];
DataTxt = [DataHeader; DataValue];

MaxCharPerCol = max(cellfun('length', DataTxt), [], 1);
MaxCharPerCol = repmat(MaxCharPerCol, size(DataTxt, 1), 1);
for j = 1:numel(DataTxt)
    DataTxt{j} = sprintf('%s%s', DataTxt{j}, repelem(' ', 1, MaxCharPerCol(j) - numel(DataTxt{j}))); 
end

LegendTxt = cell(size(DataTxt, 1), 1);
for j = 1:numel(LegendTxt)
    LegendTxt{j} = sprintf('%s  ', DataTxt{j, :});
end

X = Cutoff + 1;
Y = 0.5;
text(O.Ax.CDF, X, Y, LegendTxt, 'FontName', 'Courier New', 'FontSize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');


%Final Adjustements
setPlotTickDecimal(O.Ax.CDF, 0, 1);
title(O.Ax.CDF, [O.DataName 'CDF']);
resizeSubplots(O.Ax.CDF);

