function Gx = plotDframeUsage(Dcombo,Dname)
figure
Bx = bar(Dcombo,'stack');
Bx(1).FaceColor = [1 0 0];
Bx(2).FaceColor = [0 1 0];
Bx(3).FaceColor = [0 0 1];

Xtick = [1:size(Dcombo,1)];
Xlabel = Dname;
set(gca,'XTick',Xtick,'XTickLabel',Xlabel,'XTickLabelRotation',90)

Lx = legend({'RF 1'; 'RF 2'; 'RF 3'},'location','north');