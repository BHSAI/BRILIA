%testPlotTools tests the plotting tool functions used here.

function testPlotTools()
%Create some variables
X = 1:3:10;
Y1 = X.^2;
Y2 = X / 2;
Y3 = -X;
Gx = figure();
subplot(2,2,1) 
plot(X, Y1, 'r');
subplot(2,2,2) 
plot(X, Y2, 'g');
subplot(2,2,3) 
plot(X, Y3, 'b');
XLim = get(gca, 'XLim');
YLim = get(gca, 'YLim');
Pos = get(gca, 'Position');
AxSpecial = axes('xlim', XLim, 'ylim', YLim, 'color', 'none', 'YAxisLocation', 'right', 'XAxisLocation', 'top', 'Position', Pos);

fprintf('%s: resizeFigure(Gx, 5, 5)\n', mfilename);
resizeFigure(Gx, 5, 5);
pause(1);

fprintf('%s: setPlotTickDecimal(Gx, 2, 3)\n', mfilename);
setPlotTickDecimal(Gx, 2, 3);
pause(1);

fprintf('%s: setAxes(Gx, ''FontName'', ''Arial'', ''FontSize'', 10, ''TitleString'', ''Test'', ''TitleFontSize'', 14)\n', mfilename);
setAxes(Gx, 'FontName', 'Arial', 'FontSize', 10, 'TitleString', 'Test', 'TitleFontSize', 14);
setAxes(AxSpecial, 'XTickLabel', '', 'TitleString', '')
setPlotTickDecimal(AxSpecial, -1, 1)
pause(1);

fprintf('%s: resizeSubplots(Gx)\n', mfilename);
resizeSubplots(Gx);
pause(1);

fprintf('%s: resizeSubplots(Gx, ''VertSpacer'', 0.02, ''HorzSpacer'', 0.05)\n', mfilename);
resizeSubplots(Gx, 'VertSpacer', 0.02, 'HorzSpacer', 0.05);
pause(1);

fprintf('%s: resizeSubplots(Gx, ''FigSpacer'', 0.02)\n', mfilename);
resizeSubplots(Gx, 'FigSpacer', 0.02);
pause(1);

fprintf('%s: resizeSubplots(Gx, ''VertSpacer'', 0.02, ''HorzSpacer'', 0.05, ''FigSpacer'', 0.04)\n', mfilename);
resizeSubplots(Gx, 'VertSpacer', 0.02, 'HorzSpacer', 0.05, 'FigSpacer', 0.04);
pause(1);

fprintf('%s: labelSubplots(gcf, ''Direction'', ''Horizontal'')\n', mfilename);
labelSubplots(gcf, 'Direction', 'Horizontal')
pause(1);

fprintf('%s: labelSubplots(gcf, ''Direction'', ''Vertical'')\n', mfilename);
labelSubplots(gcf, 'Direction', 'Vertical')
pause(1);

fprintf('%s: labelSubplots(gcf, ''delete'')\n', mfilename);
labelSubplots(gcf, 'delete')
pause(1);

fprintf('%s: labelSubplots(gcf, ''Direction'', ''Horizontal'', ''FontSize'', 12))\n', mfilename);
labelSubplots(gcf, 'Direction', 'Horizontal', 'FontSize', 12)
pause(1);

fprintf('%s: labelSubplots(gcf, ''shift'', 0.01, -0.01)\n', mfilename);
labelSubplots(gcf, 'shift', 0.01, -0.01)
pause(1);

invertFigColor(gcf)
