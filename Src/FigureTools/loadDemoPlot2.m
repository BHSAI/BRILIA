function Gx = loadDemoPlot2

Gx = figure;
LineStyle = {'-', ':', '--'};
LineWidth = 1.5;

X = 1:100;
Lx = gobjects(1, 27);
Ax = gobjects(1, 9);
for f = 1:9
    Ax(f) = subplot(3, 3, f);
    for k = 1:3
        N = 3*(f-1)+k;
        Lx(N) = plot(X, 50*sin(X/10+10*rand(1)), 'LineWidth', LineWidth, 'LineStyle', LineStyle{k}, 'Color', 'k');
        if k == 1; hold on; end
        ylim([-50 50]);
        xlim([0 100]);
    end
    set(Ax(f), 'box', 'on', 'LineWidth', 1.5, 'FontName', 'Arial', 'FontSize', 10);
    
    hold off;
    if f<=3
        title('Title Here', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
    end
    if ismember(f, 1:3:9)
        ylabel({'Y label' '(deg)'}, 'FontName', 'Arial', 'FontSize', 12);
    end
    if f >= 7
        xlabel('X label', 'FontName', 'Arial', 'FontSize', 12);
    end
end

setAxes(Ax, 'XColor', [0 0 0], 'YColor', [0 0 0], 'ZColor', [0 0 0])
setAxes(Ax, 'FontSize', 12)
setAxes(Ax, 'TickLength', [.02 0.02], 'TickDir', 'both')
resizeFigure(Gx, 6.5, 4)
resizeSubplots_New(Gx, 'ScaleVertical', 'y', 'ScaleHorizontal', 'y')

%Removing X-Y labels in the interier spaces
for f = 1:length(Ax)
    if ~ismember(f, 1:3:9)
        set(Ax(f), 'YTickLabel', '')
    end
    if ~ (f >= 7)
        set(Ax(f), 'XTickLabel', '');
    end
end
