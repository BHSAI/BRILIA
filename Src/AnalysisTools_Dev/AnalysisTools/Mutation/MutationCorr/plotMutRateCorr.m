%Takes the Vmat, Dmat, Jmat from findMutationFreq, and then plots the
%dinucleotide probability for the V and the DJ segments as 12 points. A
%linear line will suggest propper annotation results. 

%Option = 'NormCol' will take the column of the 4x4 mutation matrix, and
%normalize such that each column adds to 1.
%Option = 'NormAll' will take the entire 4x4 mutation matrix, and
%normalize such that everything adds to 1.

function varargout = plotMutRateCorr(Vmat,Dmat,Jmat,varargin)
P = inputParser;
addOptional(P, 'Option', 'normcol', @(x) ischar(x) & ismember(lower(x), 'normcol', 'normall'));
addParameter(P, 'Legend', 'on', @ischar);
addParameter(P, 'Xtitle', 'on', @ischar);
addParameter(P, 'Ytitle', 'on', @ischar);
addParameter(P, 'Xlabel', 'on', @ischar);
addParameter(P, 'Ylabel', 'on', @ischar);
addParameter(P, 'FitFontSize', 18, @isnumeric);
addParameter(P, 'MarkerSize', 250, @isnumeric);

parse(P, varargin{:});
Option = P.Results.Option;
AddLegend = P.Results.Legend;
AddXtitle = P.Results.Xtitle;
AddYtitle = P.Results.Ytitle;
AddXlabel = P.Results.Xlabel;
AddYlabel = P.Results.Ylabel;
FitFontSize = P.Results.FitFontSize;
MarkerSize = P.Results.MarkerSize;

%Normalize the mutations
switch Option
  case 'normcol'
    %Normalize with respect to germline A, C, G, T
    Vmat = Vmat./repmat(sum(Vmat, 1), 4, 1);
    DJmat = Dmat+Jmat;
    DJmat = DJmat./repmat(sum(DJmat, 1), 4, 1);
  case 'normall'
    %Normalize with respect to ALL mutation pairs
    Vmat = Vmat/sum(Vmat(:));
    DJmat = Dmat+Jmat;
    DJmat = DJmat/sum(DJmat(:));
end

%Linearize and get rid of the 0 diagonal ones
Vx = Vmat(:);
Vx(1:5:end) = [];
DJy = DJmat(:);
DJy(1:5:end) = [];

%Get rid of the T's
Vx2 = Vx(1:end);
DJy2 = DJy(1:end);

%Plot the 12 dots according to the marker styles
MarkerStyles = {'d' 'o' 's' '^'}; %Germline shape
MarkerColors = {[0 0.6 0]; [1 0 0]; [.1 .1 1]; [.3 .3 .3]}; %Mutant color 
Germ2MutIdx = [...
        2 1;
        3 1;
        4 1;
        1 2;
        3 2;
        4 2;
        1 3;
        2 3;
        4 3;
        1 4;
        2 4;
        3 4]; %Left col is mutant, right side is germline.

%Set the figure size
GX = figure;
set(GX, 'units', 'pixel', 'position', [50 50 600 600]);

%Setup the plot axes properties
AX = gca;
set(AX, 'XLim', [0 1], 'YLim', [0 1]);
set(AX, 'XTick', [0:0.2:1], 'YTick', [0:0.2:1]);
set(AX, 'Box', 'on', 'TickDir', 'out', 'XGrid', 'on', 'Ygrid', 'on', 'FontSize', 16) %Fix grid and other properties
set(AX, 'Position', [0.13 0.17 0.8 0.8]);

if strcmpi(AddXtitle, 'on')
  HX = xlabel(AX, 'Mut. Freq. of V');
  set(HX, 'FontSize', 24, 'FontName', 'Arial')
end

if strcmpi(AddYtitle, 'on')
  HY = ylabel(AX, 'Mut. Freq. of DJ');
  set(HY, 'FontSize', 24, 'FontName', 'Arial')
end

%Fix the label units
if strcmpi(AddXlabel, 'on')
  Xticks = get(AX, 'XTickLabels');
  for k = 1:length(Xticks)
    Xticks{k} = sprintf('%0.1f', convStr2NumMEX(Xticks{k}));
  end
  set(AX, 'XTickLabels', Xticks);
else
  set(AX, 'XTickLabels', []);
end

if strcmpi(AddYlabel, 'on')
  Yticks = get(AX, 'YTickLabels');
  for k = 1:length(Yticks)
    Yticks{k} = sprintf('%0.1f', convStr2NumMEX(Yticks{k}));
  end
  set(AX, 'YTickLabels', Yticks);
else
  set(AX, 'YTickLabels', []);
end

resizeSubplots(gcf);

%Extend the plots to fill space
set(AX, 'OuterPosition', [0 0 1 1]);
TightInsets = get(AX, 'TightInset');
Positions = [TightInsets(1) TightInsets(2) 1-TightInsets(3)-TightInsets(1) 1-TightInsets(4)-TightInsets(2)];
set(AX, 'Position', Positions)

%Plot the data
for j = 12:-1:1
  hold(AX, 'on')
  LX = scatter(Vx(j), DJy(j), 'fill');
  hold(AX, 'off')
  Mstyle = MarkerStyles{Germ2MutIdx(j, 1)};
  Mcolor = MarkerColors{Germ2MutIdx(j, 2)};
  set(LX, 'Marker', Mstyle, 'MarkerFaceColor', Mcolor, 'SizeData', MarkerSize, 'MarkerEdgeColor', [1 1 1], 'LineWidth', 0.3)
end

%Calculate Rsq
x = Vx2;
y = DJy2;
[P, ~] = polyfit(x, y, 1);
Slope = P(1); %Difference in slope
xval = [0:0.1:1];
yval = polyval(P, xval);
hold(AX, 'on')
plot(AX, [0 1], [0 1], '--', 'color', [0.6 0.6 0.6])
plot(AX, xval, yval, '-k')

%Calculate the R pearson correlation
R = corr(x, y);
text(0.95, 0.2, sprintf('R_{corr} = %0.2f', R), 'FontName', 'Arial', 'FontSize', FitFontSize, 'HorizontalAlignment', 'right');
text(0.95, 0.1, sprintf('Slope = %0.2f', Slope), 'FontName', 'Arial', 'FontSize', FitFontSize, 'HorizontalAlignment', 'right');

%Draw the legend
if strcmpi(AddLegend, 'on')
  LegXcoor1 = [0.1 0.1 0.1 0.1];
  LegXcoor2 = LegXcoor1 + 0.15;
  LegYcoor = [0.9 0.83 0.76 0.69];
  LegText1 = {'A_{0}', 'C_{0}', 'G_{0}', 'T_{0}'};
  LegText2 = {'A_{1}', 'C_{1}', 'G_{1}', 'T_{1}'};
  hold(AX, 'on')
  for k = 1:4
    text(LegXcoor1(k)+0.03, LegYcoor(k), LegText1{k}, 'FontName', 'Arial', 'FontSize', 16, 'FontWeight', 'bold', 'VerticalAlignment', 'middle', 'Color', MarkerColors{k});
    scatter(AX, LegXcoor2(k), LegYcoor(k), 'SizeData', 250, 'Marker', MarkerStyles{k}, 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    text(LegXcoor2(k)+0.03, LegYcoor(k), LegText2{k}, 'FontName', 'Arial', 'FontSize', 16, 'VerticalAlignment', 'middle');
  end
  hold(AX, 'off')
end

if nargout >= 1
  varargout{1} = GX;
  if nargout >= 2
    varargout{2} = AX;
  end
end