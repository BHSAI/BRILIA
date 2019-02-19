%makeLogoFont generates the letter images needed for creating weblogos.
%This just needs to be run once. To change font colors, need modify
%getAATable.m
function makeLogoFont

[AAProp, AAClr] = getAATable();
AALetter = cell2mat(AAProp(:,2));

%Make sure to add the Font paths
SavePath = fullfile(fileparts(mfilename('fullpath')), 'Font');

%Generate Figure
Gx = figure('Visible', 'off');
Ax = axes(Gx);

Gx.Position = [100 100 350 300];
Gx.Units = 'points';
Gx.PaperUnits = Gx.Units;
Gx.PaperPosition = Gx.Position;

%Generate Axes
Ax.Units = 'normalized';
Ax.Position = [0 0 1 1];
Ax.XLim = [0 1];
Ax.YLim = [0 1];
axis(Ax, 'off')

%Generate fonts
Tx = text(Ax, 0.5, 0.5, '', 'FontSize', 330, 'Color', [0 0 0], 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
for j = 1:length(AALetter)
    Tx.String = AALetter(j);
    Tx.Color = AAClr(j, :);
    savePlot(Gx, 'SaveAs', fullfile(SavePath, [AALetter(j) '.png']), 'DPI', 600);    
end
close(Gx)