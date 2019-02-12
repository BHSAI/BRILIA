%scaleDistDemo demonstrates how rescaleDist works better than normalization
%by total frequency.

function rescaleDistDemo
A = [40 100 30 50  300]';
B = [40 100 30 300 300]';

Aobs = A/2;
Bobs = B*2;

Gx = figure;
resizeFigure(Gx, 4, 6)

subplot(3, 1, 1)
Bx1 = bar([Aobs, Bobs]);
YLim = get(gca, 'YLim');
text(3, YLim(2), 'Observed Freq.', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'top', 'FontName', 'Arial', 'FontSize', 12);

subplot(3, 1, 2)
Anorm = Aobs / sum(Aobs);
Bnorm = Bobs / sum(Bobs);
Bx2 = bar([Anorm, Bnorm]);
YLim = get(gca, 'YLim');
text(3, YLim(2), 'Normalized Freq.', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'top', 'FontName', 'Arial', 'FontSize', 12);

subplot(3, 1, 3)
Brs = rescaleDist(Aobs, Bobs);
Bx3 = bar([Aobs, Brs]);
YLim = get(gca, 'YLim');
text(3, YLim(2), 'Rescaled Freq.', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'top', 'FontName', 'Arial', 'FontSize', 12);

resizeSubplots;
savePlot(gcf, 'SaveAs', '');


