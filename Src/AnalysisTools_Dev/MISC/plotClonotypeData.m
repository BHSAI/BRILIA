function plotClonotypeData(ClonoData, Group)

[~, Group] = getFileList('all');

ColorMap = getStdColor;

Fields = fieldnames(ClonoData);

for p = 1:length(Fields)
    if iscell(ClonoData(1).(Fields{p})); continue; end
    figure;
    XLimVal = zeros(length(ClonoData), 1);
    for f = 1:length(ClonoData)
        Px = cdfplot(ClonoData(f).(Fields{p}));
        Px.Color = ColorMap(Group(f), :);
        XLimVal(f) = Px.XData(find(Px.YData > 0.95, 1));
        if f == 1; hold on; end
    end
    hold off
    title(Fields{p});
    XLim = [0 3*max(XLimVal)];
    YLim = get(gca, 'YLim');
    xlim(XLim);
    text(gca, XLim(2), YLim(1)+ 3*0.05*YLim(2), 'Naive  ', 'Color', ColorMap(1, :), 'HorizontalAlignment', 'right', 'FontSize', 14, 'FontName', 'Arial');
    text(gca, XLim(2), YLim(1)+ 2*0.05*YLim(2), 'VLP  ', 'Color', ColorMap(2, :), 'HorizontalAlignment', 'right', 'FontSize', 14, 'FontName', 'Arial');
    text(gca, XLim(2), YLim(1)+ 1*0.05*YLim(2), 'VLP+ADJ  ', 'Color', ColorMap(3, :), 'HorizontalAlignment', 'right', 'FontSize', 14, 'FontName', 'Arial');    
    setPlotTickDecimal(gca, 0, 1);
    resizeSubplots(gcf);
    savePlot(gcf, 'SaveAs', sprintf('%s.png', Fields{p}));
end
