%makeLogoFont generates the letter images needed for creating weblogos.
%This just needs to be run once. To change font colors, need modify
%getAATable.m
function makeLogoFont
%Determine the Font image folder
SavePath = fullfile(fileparts(mfilename('fullpath')), 'Font');

%Get the amino acid properties and color
[AAProp, AAClr] = getAATable();
AALetter = cell2mat(AAProp(:,2));

%Generate Figure
Gx = figure('Visible', 'on');
Ax = axes(Gx);

Gx.Position = [100 100 340 300];
Gx.PaperUnits = Gx.Units;
Gx.PaperPosition = Gx.Position;

%Generate Axes
Ax.Units = 'normalize';
Ax.XLim = [0 1];
Ax.YLim = [0 1];
Ax.Position = [0 0 1 1];
axis(Ax, 'on')

%Generate fonts
%Determine the smalles font height for rescaling purposes

Gx = figure;
Tx = annotation('textbox',[0 0 1 1], 'String', 'C', 'FontSize', 300, 'FontName', 'Lucida Console', 'Color', [0 0 0], 'FitBoxToText', 'on', 'LineStyle', 'none');
Gx.Units = 'inches';
Tx.Units = 'inches';

Gx.Units = 'points';
Ax = axes(Gx);
Ax.Units = 'points';
Tx = text(Ax, 0, 1, 'C', 'FontSize', 300, 'Color', [0 0 0], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'cap', 'FontName', 'Lucida Console');


Tx.Units = 'inches';

W = Tx.Extent(1) + Tx.Extent(3);
H = Tx.Extent(2) + Tx.Extent(4);
Gx.Position(3:4) = [W H];
Ax.Position = [0 0 W H];
Tx.Position(1:2) = 0;


Ax.Position(3:4) = Gx.Position(3:4);
Ax.Position(1:2) = 0;
Tx = text(Ax, 0.5, 0.5, 'C', 'FontSize', 300, 'Color', [0 0 0], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontName', 'Lucida Console');
Gx.Position(3) = Tx.Extent(1) + Tx.Extent(3);
Gx.Position(4) = Tx.Extent(2) + Tx.Extent(4);
Ax.Position(3) = Gx.Position(3)
Ax.Position(4) = Gx.Position(4);
Gx.OuterPosition(3:4) = Gx.Position(3:4);

Tx.Units = 'inches';
Ax.Units = Gx.Units;

Ax.Position(3:4) = Gx.Position(3:4);
Ax.Position(1:2) = 0;


Tx.Extent
Gx.Position(3) = Tx.Extent(2) + Tx.Extent(4)

Ax.Units = 'points';
Ax.Position(3) = Tx.Extent(3) - Tx.Extent(1)
Ax.Position(4) 
Gx.Position(3:4) = Tx.Extent(3:4)

for j = 1:numel(AALetter)
    Tx.String = AALetter(j);
    Tx.Color = AAClr(j, :);
    savePlot(Gx, 'SaveAs', fullfile(SavePath, [AALetter(j) '.png']), 'DPI', 600);    
end
close(Gx)