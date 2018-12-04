%plotTree will plot the lineage trees for each cluster > 1 member. All
%plots will be saved to a folder named "plotTree", which is in the same
%directory where the BRILIA csv file is, or where it is specified using the
%('SaveAs', [filename]) parameter pair.
%
%  plotTree
%
%  plotTree(FileName)
%
%  plotTree(FileName, P)
%
%  plotTree(VDJdata, VDJheader, Param, Value, ...)
%
% [ImageGrpNums, ImageNames, TreeSearchTable] = plotTree(...)
%
%  INPUT
%    VDJdata: BRILIA output table
%    VDJheader: BRILIA header cell
%    FileName: file name of the BRILIA csv output
%    P: a structure containing P.(Param) = Value information, obtained from
%      P = plotTree('getinput');
%
%      Parmaeter Name  Value       Description
%      --------------- ---------   ----------------------------------------
%      DistanceUnit    'ham'       Set x-axis to hamming distance.
%                      'hamperc'   Set x-axis to hamming distance %.
%      DotMaxSize      N >= 100    Max dot area size, in (pixel^2). Default
%                                    is 300.
%      DotMinSize      N >= 1      Min dot area size, in (pixel^2). Default
%                                    is 1.
%      DotScalor       N >= 1      Multiply template count by this to get
%                                    DotSize per each sequence. Will not
%                                    exceed DotMax. Default is 30.
%      DotColorMap     Mx3 matrix  A RGB colormap matrix used to decide how
%                                    to color each tree node. Default is a 
%                                    jet colormap.
%      GrpMinSize      N >= 0      Minimum tree size to draw
%      GrpMaxSize      N >= 1      Maximum tree size to draw
%      Legend          'y', 'n'    Will draw a color-coded CDR3 legend.
%                                    Default is y.
%      LegendFontSize  10          Font size of the Legend
%      Xmax            N > 0       X axis maximum value.
%      Yincr           N > 0       Vertical spacing between nodes. Default
%                                    is 0.125".
%      FigWidth        3.3         Width of the whole figure, inch
%      FigMaxHeight    5           Maximum height of the whole fig, inch
%      SaveAs          filepath    Will save to this folder, using default
%                                    file names as Tree.GrpN.png
%                      filename    Will save to this file's folder(if any),
%                                    and file names are SaveAs.GrpN.png
%                                    IF no path is specified in SaveAs and
%                                    input file is known, will save to that
%                                    directory plotTree subfolder.
%      ShowStatus      'y' 'n'     If y, will show how many trees are left
%                                    for saving.
%
%  OUTPUT
%   ImageGrpNums: group numbers that have tree images
%   ImageNames: image file names in the order of the ImageGrpNums 
%   TreeSearchTable: simplified table to look at tree image file and group
%     information.
%
%  NOTE
%    To use a custom color scheme for each tree node, assemble the Mx3 RGB
%    matrix separately and then use that as a DotColor input.
% 
%  EXAMPLE 
%    plotTree('briliaoutput.csv', 'SaveAs', '/newdir/plotTree/tree.png')

function varargout = plotTree(varargin)
varargout = cell(1, nargout);

P = inputParser;
addParameter(P, 'GrpNum',           [],         @(x) isnumeric(x) && x >= 1);
addParameter(P, 'GrpMinSize',       1,          @(x) isnumeric(x) && x >= 1);
addParameter(P, 'GrpMaxSize',       Inf,        @(x) isnumeric(x) && x >= 1);
addParameter(P, 'DistanceUnit',     'hamperc',  @(x) ischar(x) && any(ismember(lower(x), {'shm', 'ham', 'shmperc', 'hamperc'})));
addParameter(P, 'DotMaxSize',       300,        @(x) isnumeric(x) && x >= 1);
addParameter(P, 'DotMinSize',       1,          @(x) isnumeric(x) && x >= 1);
addParameter(P, 'DotScalor',        5,          @(x) isnumeric(x) && x >= 1);
addParameter(P, 'DotColorMap',      [],         @(x) isempty(x) || (isnumeric(x) && size(x, 2) == 3));
addParameter(P, 'Legend',           'y',        @(x) ischar(x) && any(ismember(lower(x), {'y', 'n'})));
addParameter(P, 'LegendFontSize',   10,         @(x) isnumeric(x) && x >= 1);
addParameter(P, 'Xmax',             [],         @(x) isempty(x) || (isnumeric(x) && x > 10));
addParameter(P, 'Yincr',            0.125,      @(x) isnumeric(x) && x > 0);
addParameter(P, 'FigMaxHeight',     5,          @(x) isnumeric(x) && x > 1);
addParameter(P, 'FigWidth',         3.3,        @(x) isnumeric(x) && x > 1);
addParameter(P, 'FigSpacer',        0.01,       @(x) isnumeric(x) && x >= 0)
addParameter(P, 'Visible',          'off',      @(x) ischar(x) && ismember(lower(x), {'on', 'off'}));
addParameter(P, 'SaveAs',           '',         @(x) ischar(x) || isempty(x));
addParameter(P, 'SaveDir',          '',         @(x) ischar(x) || isempty(x));
addParameter(P, 'SaveSubDir',       'Tree',     @(x) ischar(x) || isempty(x));
addParameter(P, 'StatusHandle',     [],         @(x) isempty(x) || ishandle(x) || strcmpi(class(x), 'matlab.ui.control.UIControl'));
addParameter(P, 'ShowOnly',         false,      @islogical);

if ~(~isempty(varargin) && ischar(varargin{1}) && ismember(lower(varargin{1}), {'getinput', 'getinputs', 'getvarargin'}))
    [VDJdata, VDJheader, ~, FilePath, varargin] = getPlotVDJdata(varargin{:});
    if isempty(VDJdata); return; end
end

[Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInput(P, varargin{:});
if ReturnThis
   varargout = {Ps, Pu, ExpPs, ExpPu};
   return
end

GetGrpNum = Ps.GrpNum;
GrpMinSize = Ps.GrpMinSize;
GrpMaxSize = Ps.GrpMaxSize;
DistanceUnit = Ps.DistanceUnit;
DotMaxSize = Ps.DotMaxSize;
DotMinSize = Ps.DotMinSize;
DotScalor = Ps.DotScalor;
DotColorMap = Ps.DotColorMap;
Legend = Ps.Legend;
LegendFontSize = Ps.LegendFontSize;
Xmax = Ps.Xmax;
Yincr = Ps.Yincr;
FigMaxHeight = Ps.FigMaxHeight;
FigWidth = Ps.FigWidth;
FigSpacer = Ps.FigSpacer;
Visible = Ps.Visible;
SaveAs = Ps.SaveAs;
SaveDir = Ps.SaveDir;
SaveSubDir = Ps.SaveSubDir;
StatusHandle = Ps.StatusHandle;
ShowOnly = Ps.ShowOnly;

%Check if Min/Max values are valid
if DotMinSize > DotMaxSize
    warning('%s: DotMinSize cannot be > DotMaxSize. Setting DotMinSize = DotMaxSize', mfilename)
    DotMinSize = DotMaxSize;
end

if GrpMinSize > GrpMaxSize
    warning('%s: GrpMinSize cannot be > GrpMaxSize. Setting GrpMinSize = GrpMaxSize', mfilename)
    GrpMinSize = GrpMaxSize;
end

%Determine where to save
if ~ShowOnly
    if isempty(SaveAs)
        SaveAs = prepSaveTarget('SaveDir', FilePath, 'SaveSubDir', SaveSubDir, 'SaveAs', 'Tree.png', 'MakeSaveDir', 'y');
    else
        SaveAs = prepSaveTarget('SaveDir', SaveDir, 'SaveSubDir', SaveSubDir, 'SaveAs', SaveAs, 'MakeSaveDir', 'y');
    end
else
    SaveAs = '';
    Visible = 'on';
end

%--------------------------------------------------------------------------

%Select relevant groups and splice VDJdata for parfor
G = getGrpIdx(VDJdata, VDJheader, GetGrpNum);
KeepLoc = [G.Size] >= GrpMinSize & [G.Size] <= GrpMaxSize;
G = G(KeepLoc);
Map = getVDJmapper(VDJheader);
if isempty(G); return; end
VDJdata = spliceData(VDJdata([G.Idx], :), VDJheader);

%For storing the actual file names
ImgNames = cell(length(G), 1);
ImgGrpNums = zeros(length(G), 1);

%For tracking progress
PT = ProgressTracker(length(G), [], 'Trees drawn', StatusHandle);
DQ = PT.DataQueue;
showStatus(sprintf('Drawing %d trees.', length(G)), StatusHandle)

if ShowOnly %Do not use parfor to preserve figure handles
    for y = 1:length(G)
        try
            ExpPuT = ExpPu;
            VDJdata{y} = padtrimSeqGroup(VDJdata{y}, Map, 'grpnum', 'trim', 'Seq'); 
            plotSingleTree(VDJdata{y}, VDJheader, ExpPuT, DistanceUnit, DotMaxSize, DotMinSize, DotScalor, DotColorMap, Legend, LegendFontSize, Xmax, Yincr, FigMaxHeight, FigWidth, FigSpacer, Visible);
        catch ME
            throw(ME)
        end
    end
    return
else
    parfor y = 1:length(G)
        try
            ExpPuT = ExpPu;
            VDJdata{y} = padtrimSeqGroup(VDJdata{y}, Map, 'grpnum', 'trim', 'Seq'); 
            Gx = plotSingleTree(VDJdata{y}, VDJheader, ExpPuT, DistanceUnit, DotMaxSize, DotMinSize, DotScalor, DotColorMap, Legend, LegendFontSize, Xmax, Yincr, FigMaxHeight, FigWidth, FigSpacer, Visible);

            SavePre = sprintf('.Grp%d', G(y).GrpNum);
            SuggestedSaveName = prepSaveTarget('SaveAs', SaveAs, 'SavePrefix', SavePre);
            ImgNames{y} = savePlot(Gx, 'SaveAs', SuggestedSaveName, ExpPuT{:});
            ImgGrpNums(y) = G(y).GrpNum;
            close(Gx);
        catch ME
            throw(ME)
        end
        send(DQ, y); %To show progress
    end
end

%Create the tree search table
VDJdata = joinData(VDJdata, VDJheader, 'stable');
TreeSearchTable = getTreeSearchTable(VDJdata, VDJheader);
[~, ~, TableIdx] = intersect(ImgGrpNums, cell2mat(TreeSearchTable(2:end, 1)));
TreeSearchTable = TreeSearchTable([1; TableIdx+1], :);
for j = 1:length(ImgNames) %Remove the NewFilePath from the TreeFiles
    [FilePath, FileName, ~] = parseFileName(ImgNames{j});
    ImgNames{j} = FileName;
end
TreeSearchTable(2:end, end) = ImgNames;
writeDlmFile(TreeSearchTable, fullfile(FilePath, 'TreeSearchTable.csv'));
varargout{1} = ImgGrpNums;
varargout{2} = ImgNames;
varargout{3} = TreeSearchTable;

function Gx = plotSingleTree(VDJdata, VDJheader, ExpPu, DistanceUnit, DotMaxSize, DotMinSize, DotScalor, DotColorMap, Legend, LegendFontSize, Xmax, Yincr, FigMaxHeight, FigWidth, FigSpacer, Visible)

[AncMapS, TreeName, CDR3Name, TemplateCount] = getTreeData(VDJdata, VDJheader);
AncMap = AncMapS.(upper(DistanceUnit));
TreeCoord = calcTreeCoord(AncMap);
if isempty(TreeCoord); return; end

MaxHorzCoord = max(TreeCoord(:, 1));
MaxVertCoord = max(TreeCoord(:, 2));

%Determine the legend and colors
if strcmpi(Legend, 'y') && ~isempty(CDR3Name)
    [CDR3legend, UnqCDR3seq] = makeTreeLegend_CDR3(CDR3Name);
    [DotColor, UnqDotColor] = mapDotColor_CDR3(CDR3Name, UnqCDR3seq, 'ColorMap', DotColorMap);
end

[Gx, Ax] = makeDefaultFigure(Visible, FigWidth, FigMaxHeight, DistanceUnit, ExpPu);    
cla(Ax); %Clear and add title now. Required early to get the axes tight inset.
set(get(Ax, 'Title'), 'String', TreeName);
drawnow('nocallbacks');

%----------------------------------------------------------------------
%Calculate the figure and axes heights

TreeHeight = Yincr * (MaxVertCoord + 1);

AxesBorder = get(Ax, 'TightInset') + ones(1, 4)*FigSpacer;
AxesBorder([1 3]) = max(AxesBorder([1 3]));
AxesBorderHeight = sum(AxesBorder([2 4]));

LegendHeight = 0;
if strcmpi(Legend, 'y')
    LegendFontHeight = (LegendFontSize + 1)/72; %Add a 1pt spacer
    LegendHeight = LegendFontHeight * length(UnqCDR3seq);
end

AxesHeight = max(LegendHeight, TreeHeight);
FigHeight = min(AxesHeight + AxesBorderHeight, FigMaxHeight);
if FigHeight == FigMaxHeight
    AxesHeight = FigMaxHeight - AxesBorderHeight;
end

Gx.Position(4) = FigHeight;
resizeSubplots(Gx, 'FigSpacer', FigSpacer);

%----------------------------------------------------------------------
%Calculate the X and Y limits.

Xlim = [0, MaxHorzCoord];
Ylim = [0, MaxVertCoord + 1];  

%Calculate the XLim required to fit the legend without overlapping
if strcmpi(Legend, 'y') && ~isempty(CDR3legend{1})
    ReqLegendTextHeight = min(AxesHeight/length(CDR3legend), (LegendFontSize + 1) / 72);         
    ReqLegendFontSize = round(ReqLegendTextHeight * 72) - 1; %Remember, 1 pt spacer rule
    if ReqLegendTextHeight < (LegendFontSize + 1) / 72
        ReqLegendFontSize = max(2, ReqLegendFontSize);
        ReqLegendTextHeight = (ReqLegendFontSize + 1) / 72;
    end

    %Figure out what the text width would be
    HorzSpacer = 0.01; %inches, adding 0.01 in right side spacer
    TempText = text(1, 1, CDR3legend{1}, 'FontName', 'Courier', 'FontSize', ReqLegendFontSize, 'Units', 'inches'); %Use courier for even spacing. Include spacing at end.
    TextExt = get(TempText, 'Extent');
    TextWidth = TextExt(3) + HorzSpacer; 
    delete(TempText)

    Xlim(2) = Xlim(1) + (MaxHorzCoord - Xlim(1)) * Ax.Position(3) / (Ax.Position(3) - TextWidth); %Remember, dots can bleed into text.
end

%----------------------------------------------------------------------
%Retain dots inside axes bounds by adjusting Ylim, Xlim, or DotScalor

EdgePerc = 0.05; %Amount of left and right edge to preserve for X limit
DotSizes = TemplateCount * DotScalor;
DotRadii = (DotSizes / pi).^0.5 / 72;
InchPerX = (1 - EdgePerc*2) * Ax.Position(3) / diff(Xlim); %Remember, 5% of fig width are used for rescaling xlim

%Rescale xlim(2) if dots overlap with legend or right border
DotRights = InchPerX * (TreeCoord(:, 1) - Xlim(1)) + EdgePerc * Ax.Position(3) + DotRadii;
if strcmpi(Legend, 'y') && ~isempty(CDR3legend)
    RightEdge = Ax.Position(3) - TextWidth;
else
    RightEdge = Ax.Position(3);
end
ExpandRight = max(DotRights) - RightEdge;
if ExpandRight > 0
    Xlim(2) = ((1 - EdgePerc) * Ax.Position(3) * Xlim(2) - ExpandRight * Xlim(1)) / ((1 - EdgePerc * 2) * Ax.Position(3) - ExpandRight);
    InchPerX = (1 - EdgePerc*2) * Ax.Position(3) / diff(Xlim);
end

%Checking how much to rescale dot due to overflow in left side
DotLefts = InchPerX * (TreeCoord(:, 1) - Xlim(1)) + EdgePerc * Ax.Position(3) - DotRadii;
ExpandLeft = -min(DotLefts);
if ExpandLeft > 0 %Go straight for rescaling dot sizes
    LeftestDot = find(DotLefts == min(DotLefts));
    LeftestDotRad = DotRadii(LeftestDot(1));
    DotRadii = DotRadii * (LeftestDotRad - ExpandLeft) / LeftestDotRad;
end

InchPerY = Ax.Position(4) / diff(Ylim);

%Rescale the axes height, if possible if dots go over top/bot
AvailExpHeight = FigMaxHeight - (AxesBorderHeight + Ax.Position(4));
DotHighs = InchPerY * (TreeCoord(:, 2) - Ylim(1)) + DotRadii;
ExpandTop = max(DotHighs) - Ax.Position(4);
if ExpandTop > 0 
    HighestDotIdx = find(DotHighs == max(DotHighs), 1);
    HighestDotRad = DotRadii(HighestDotIdx);
    if ExpandTop < AvailExpHeight
        Gx.Position(4) = Gx.Position(4) + ExpandTop;
        Ax.Position(4) = Ax.Position(4) + ExpandTop;
        Ylim(2) = Ylim(2) + ExpandTop / InchPerY;
        AvailExpHeight = AvailExpHeight - ExpandTop;
    else
        DotRadii = DotRadii * (HighestDotRad - ExpandTop) / HighestDotRad; 
    end
end

DotLows = InchPerY * (TreeCoord(:, 2) - Ylim(1)) - DotRadii;
ExpandBot = -min(DotLows);
if ExpandBot > 0 
    LowestDotIdx = find(DotLows == min(DotLows), 1);
    LowestDotRad = DotRadii(LowestDotIdx);
    if ExpandBot < AvailExpHeight
        Ax.Position(4) = Ax.Position(4) + ExpandBot;
        Ylim(1) = Ylim(1) - ExpandBot / InchPerY;
    else
        DotRadii = DotRadii * (LowestDotRad - ExpandBot) / LowestDotRad;
    end
end

DotSizes = round(pi * (DotRadii * 72).^2); % Final dot size to use
DotSizes(DotSizes == 0) = 1;
if TemplateCount(1) == 0 %Ensure inferred root DotSize is 1. Inferred Root has a 0 template count.
    DotSizes(1) = 1;
end

%Rescale dots based on custom min/max sizes. WARNING: This can cause
%dots to go out of bounds of the plots.
RealMinDotSize = min(DotSizes(DotSizes > 0));
if RealMinDotSize < DotMinSize %Ensure lower DotSize >= DotMinSize
    DotSizes = DotSizes + (DotMinSize - RealMinDotSize);
end
RealMaxDotSize = max(DotSizes(DotSizes > 0));
if RealMaxDotSize > DotMaxSize %Ensure upper DotSize <= DotMaxSize
    DotScaleFactor = (DotMaxSize - DotMinSize) / (RealMaxDotSize - DotMinSize);
    DotSizes = (DotSizes - DotMinSize) * DotScaleFactor + DotMinSize;
end

%Center vertically
InchPerY = Ax.Position(4) / diff(Ylim);
DotHighs = InchPerY * (TreeCoord(:, 2) - Ylim(1)) + DotRadii;
DotLows = InchPerY * (TreeCoord(:, 2) - Ylim(1)) - DotRadii;
TopGap = Ax.Position(4) - max(DotHighs);
BotGap = min(DotLows);
if TopGap > BotGap
    Yshift = -(((TopGap + BotGap) / 2) - BotGap) / InchPerY; %Negative means shift up
else
    Yshift = +(((TopGap + BotGap) / 2) - TopGap) / InchPerY; %Positive means shift down
end
Ylim(1) = Ylim(1) + Yshift;

%Adjust Xlim for aesthetically good locations
if ~isempty(Xmax)
    Xlim(2) = Xmax;
end
if Xlim(2) < 10
    Xlim(2) = 10;
end
Xlim(2) = ceil(Xlim(2) / 5) * 5; %Ensure ends at multiple of 5
Xlim = [(-EdgePerc*Xlim(2))  ((1+EdgePerc)*Xlim(2))]; %Ensure edges are constant 5% off the edge. Prevents X labels from going over too.

%You want at most 10 values on X, and multiple of 5's if possible)
Xincr = ceil(ceil(Xlim(2)/10)/5)*5;
XTickVal = 0:Xincr:Xlim(2);
XTickLab = cell(size(XTickVal));
for w = 1:length(XTickLab)
    XTickLab{w} = num2str(XTickVal(w), '%d');
end

%Set axes prop BEFORE plotting to avoid some cutoff issues.
set(Ax, 'XLim', Xlim, ...
        'YLim', Ylim, ...
        'XTick', XTickVal, ...
        'XTickLabel', XTickLab, ...
        'XMinorTick', 'off', ...
        'YTickLabel', '');

%Draw tree lines first
hold(Ax, 'on')
for j = 1:size(TreeCoord, 1)
    ParentLoc = AncMap(j, 2);
    if ParentLoc == 0 %if it's the root, use TreeCoord value directly.
        ParentLoc = j;
    end
    X0 = TreeCoord(ParentLoc, 1);
    Y0 = TreeCoord(ParentLoc, 2);
    X1 = TreeCoord(j, 1);
    Y1 = TreeCoord(j, 2);
    plot([X0 X0 X1], [Y0 Y1 Y1], 'k', 'LineWidth', 1)
end

%Draw tree nodes and leaves next
[~, IdxT] = sort(DotSizes, 'descend'); %Sort to prevent covering up smaller dot
if strcmpi(Legend, 'y') && ~isempty(CDR3Name)
    scatter(Ax, TreeCoord(IdxT, 1), TreeCoord(IdxT, 2), DotSizes(IdxT, :), DotColor(IdxT, :), 'fill');
else
    scatter(Ax, TreeCoord(IdxT, 1), TreeCoord(IdxT, 2), DotSizes(IdxT, :), [0 0 0], 'fill');
end
hold(Ax, 'off')

%Draw the legend
if strcmpi(Legend, 'y') && ~isempty(CDR3legend)
    %XY coordinate, right-aligned. Use inch, with resect to plot botleft corner.
    PlotPosition = get(Ax, 'Position');
    XYcoor = zeros(size(CDR3legend, 1), 2);
    if ReqLegendTextHeight * length(CDR3legend) < PlotPosition(4)
        XYcoor(:, 2) = PlotPosition(4) - ReqLegendTextHeight * (1:length(CDR3legend)); %Vertical anchor points
    elseif length(CDR3legend) == 1
        XYcoor(1, 2) = PlotPosition(4) - ReqLegendTextHeight;
    else
        XYcoor(:, 2) = (PlotPosition(4) - ReqLegendTextHeight) - (PlotPosition(4) - ReqLegendTextHeight) / (length(CDR3legend) - 1) * (0:length(CDR3legend) - 1);
    end
    XYcoor(:, 1) = PlotPosition(3) - HorzSpacer; %Horizontal anchor points, based from the right border of plot. With 0.01 in spacer.

    %Add the legend text on the plots
    for j = 1:length(CDR3legend)
        TextName = strrep(CDR3legend{j}, '_', '\_'); %Ensures any underscores prevent making subscripts
        text(XYcoor(j, 1), XYcoor(j, 2), TextName, 'FontWeight', 'bold', 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'bottom', 'FontName', 'Courier', 'FontSize', ReqLegendFontSize, 'Color', UnqDotColor(j, :), 'Units', 'inch');
    end
end


function varargout = makeDefaultFigure(Visible, FigWidth, FigMaxHeight, DistanceUnit, ExpPu)
DefaultFormat = {'FontName', 'Arial', ...
                 'FontSize', 10, ...
                 'TitleFontName', 'Arial', ...
                 'TitleFontSize', 12, ...
                 'TickLength', [0.005 0.005], ...
                 'TickDir', 'both' ...
                 'YTickLabelMode', 'manual', ...
                 'YTickLabel', '',  ...
                 'XTickLabelMode', 'auto', ...
                 'Box', 'on' ...
                 'Units', 'inches'}; 
Gx = resizeFigure('FigWidth', FigWidth, 'FigHeight', FigMaxHeight, 'Visible', Visible, 'Recenter', 'n');
set(Gx, 'renderer', 'painters', 'Units', 'inches');
Ax = get(Gx, 'children');
XaxisName = strrep(upper(DistanceUnit), 'PERC', ' %');
xlabel(Ax, XaxisName);
setAxes(Ax, DefaultFormat{:});
setAxes(Ax, ExpPu{:}); %Any axes setting will override defaults
drawnow('nocallbacks');

varargout{1} = Gx;
varargout{2} = Ax;