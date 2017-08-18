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
%  plotTree(Param, Value, ...)
%
%  [ImageGrpNums, ImageNames, TreeSearchTable] = plotTree(...)
% 
%  P = plotTree('getinput')
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
%                      'shm'       Set x-axis to shm distance.
%                      'shmperc'   Set x-axis to shm distance %.
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
%      Legend          'y', 'n'     Will draw a color-coded CDR3 legend.
%                                    Default is y.
%      LegendFontSize  10          Font size of the Legend
%      Sort            'y', 'n'     Will sort the tree to cluster related
%                                    branches together.
%      Xmax            N > 0       X axis maximum value. Default 20.
%      Yincr           N > 0       Vertical spacing between nodes. Default
%                                    is 0.125".
%      FigWidth        3.3         Width of the whole figure, inch
%      FigMaxHeight    5           Maximum height of the whole fig, inch
%      SaveAs          filepath    Will save to this folder, using default
%                                    file names as Tree.GrpN.png
%                      filename    Will save to this file's folder(if any),
%                                    and file names are SaveAs.GrpN.png
%                                  IF no path is specified in SaveAs and
%                                    input file is known, will save to that
%                                    directory plotTree subfolder.
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
P = inputParser;
addParameter(P, 'GetGrpNum', [], @(x) isempty(x) || isnumeric(x));
addParameter(P, 'GetSeqNum', [], @(x) isempty(x) || isnumeric(x));
addParameter(P, 'GetSizeRange', [], @(x) isempty(x) || isnumeric(x));
addParameter(P, 'GetCDR3seq', [], @(x) isempty(x) || ischar(x) || iscell(x));
addParameter(P, 'DistanceUnit', 'hamperc', @(x) any(validatestring(lower(x), {'shm', 'ham', 'shmperc', 'hamperc'})));
addParameter(P, 'DotMaxSize', 300, @(x) isnumeric(x) && x >= 1);
addParameter(P, 'DotMinSize', 1, @(x) isnumeric(x) && x >= 1);
addParameter(P, 'DotScalor', 5, @(x) isnumeric(x) && x >= 1);
addParameter(P, 'DotColorMap', [], @(x) isempty(x) || (isnumeric(x) && size(x, 2) == 3));
addParameter(P, 'Legend', 'y', @(x) any(validatestring(lower(x), {'y', 'n'})));
addParameter(P, 'LegendFontSize', 10, @isnumeric)
addParameter(P, 'Sort', 'y', @(x) any(validatestring(lower(x), {'y', 'n'})));
addParameter(P, 'Xmax', 0, @isnumeric);
addParameter(P, 'Yincr', 0.125, @isnumeric);
addParameter(P, 'FigMaxHeight', 5, @isnumeric);
addParameter(P, 'FigWidth', 3.3, @isnumeric);
addParameter(P, 'Visible', 'off', @(x) ischar(x) && ismember(lower(x), {'on', 'off'}));
addParameter(P, 'FigSpacer', 0, @(x) isnumeric(x) && x >= 0)
addParameter(P, 'SaveAs', '', @ischar);
addParameter(P, 'SaveDir', '', @(x) ischar(x) || isempty(x));
addParameter(P, 'SaveSubDir', 'Analysis', @(x) ischar(x) || isempty(x));

if ~(~isempty(varargin) && ischar(varargin{1}) && ismember(lower(varargin{1}), {'getinput', 'getinputs', 'getvarargin'}))
    [VDJdata, VDJheader, ~, FilePath, varargin] = getPlotVDJdata(varargin{:});
    if isempty(VDJdata)
        varargout = cell(1, nargout);
        return;
    end
    [H, ~, ~] = getAllHeaderVar(VDJheader);
    GrpNum = cell2mat(VDJdata(:, H.GrpNumLoc));
    UnqGrpNum = unique(GrpNum);
    %VDJdata = filterVDJdata(VDJdata, VDJheader, varargin{:});
end

[Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInput(P, varargin{:});
if ReturnThis
   varargout = {Ps, Pu, ExpPs, ExpPu};
   return;
end

%Check if P.DotMinSize <= P.DotMaxSize
if Ps.DotMinSize > Ps.DotMaxSize
    Ps.DotMinSize = Ps.DotMaxSize;
end

%==========================================================================
%Determine where to save
if isempty(Ps.SaveAs)
    FullSaveName = prepSaveTarget('SaveDir', FilePath, 'SaveSubDir', Ps.SaveSubDir, 'SaveAs', 'Tree.png', 'MakeSaveDir', 'y');
    Ps.SaveAs = FullSaveName;
else
    FullSaveName = prepSaveTarget('SaveDir', Ps.SaveDir, 'SaveSubDir', Ps.SaveSubDir, 'SaveAs', Ps.SaveAs, 'MakeSaveDir', 'y');
    Ps.SaveAs = FullSaveName;
end
ImageNames = cell(length(UnqGrpNum), 1);
ImageGrpNums = zeros(length(UnqGrpNum), 1);

%Create a default figure
DefaultFormat = {'FontName', 'Arial', 'FontSize', 10, ...
                 'TitleFontName', 'Arial', 'TitleFontSize', 12, ...
                 'TickLength', [0.005 0.005], 'TickDir', 'both' ...
                 'YTickLabelMode', 'manual', 'YTickLabel', '', 'YTick', [], ...
                 'XTickLabelMode', 'auto', ...
                 'box', 'on'};
Gx = figure('Visible', Ps.Visible);
Ax = axes(Gx);
pause(0.05); %Wait 50ms for concurrency issues
Gx = resizeFigure(Gx, 'FigWidth', Ps.FigWidth, 'FigHeight', Ps.FigMaxHeight, 'Recenter', 'y');
XaxisName = [strrep(upper(Ps.DistanceUnit), 'PERC', ' %') ' Distance'];
xlabel(Ax, XaxisName);
set(Ax, 'Units', 'inches');
setAxes(Ax, DefaultFormat{:}); 
setAxes(Ax, ExpPu{:}); %Any axes setting will override defaults

for y = 1:length(UnqGrpNum)
    if mod(y, 100) == 0
        fprintf('%s: Trees processed: %d of %d.\n', mfilename, y, length(UnqGrpNum));
    end
    
    %Isolate group data and name
    GrpIdx = GrpNum == UnqGrpNum(y);
    if sum(GrpNum) <= 1; continue; end
    Tdata = VDJdata(GrpIdx, :);
    [AncMapS, TreeName, CDR3Name, TemplateCount] = getTreeData(Tdata, VDJheader);
    AncMap = AncMapS.(upper(Ps.DistanceUnit));
    TreeCoord = calcTreeCoord(AncMap);    
    
    %Make the legend text and get default color scheme using CDR3
    [CDR3legend, UnqCDR3seq] = makeTreeLegend_CDR3(CDR3Name);
    [DotColor, UnqDotColor] = mapDotColor_CDR3(CDR3Name, UnqCDR3seq, 'ColorMap', Ps.DotColorMap);
    
    %Get dot sizes and resize dots if they exceed the min or max dot size.
    DotSize = TemplateCount*Ps.DotScalor;
    if min(DotSize(DotSize > 0)) < Ps.DotMinSize
        DotSize = DotSize + Ps.DotMinSize - min(DotSize(DotSize >0));
    end
    if max(DotSize) > Ps.DotMaxSize
        DotSize = DotSize * Ps.DotMaxSize / max(DotSize);
    end
    DotSize(TemplateCount <= 0) = 1;
    
    %----------------------------------------------------------------------
    %Drawing the lineage tree

    cla(Ax);
    delete(get(Ax, 'title'))
    title(Ax, TreeName);

    %Draw the lineage lines only
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

    %Draw tree nodes and leaves
    [~, IdxT] = sort(DotSize, 'descend'); %Sort to prevent covering up smaller dot
    scatter(Ax, TreeCoord(IdxT, 1), TreeCoord(IdxT, 2), DotSize(IdxT, :), DotColor(IdxT, :), 'fill')
    hold(Ax, 'off')

    %Reset the labels and axis lengths
    MaxX = max(TreeCoord(:, 1));
    MaxY = max(TreeCoord(:, 2));
    set(Ax, 'Ylim', [0 MaxY+1]);
    set(Ax, 'Xlim', [-2, MaxX+2]);
    
    %----------------------------------------------------------------------
    %Rescaling vertical components depending on legend

    %Determine what the figure height should be
    LegendTextHeight = (Ps.LegendFontSize + 1)/72;
    TotalLegendHeight = length(CDR3legend) * LegendTextHeight;
    if Ps.Yincr == 0
        Ps.Yincr = TextHeight;
    end
    AxesHeight = (MaxY+1) * Ps.Yincr; %Height of all tree
    AxesHeight = max(AxesHeight, TotalLegendHeight);
    
    %IF FigHeight > FigMaxHeight, recalc height of axes, figure, and legend 
    AxesTights = get(Ax, 'TightInset');
    VertTights = sum(AxesTights([2 4]));
    FigHeight = AxesHeight + VertTights;
    if FigHeight > Ps.FigMaxHeight
        ScaleFactor = (Ps.FigMaxHeight - VertTights) / (FigHeight - VertTights);
        AxesHeight = AxesHeight * ScaleFactor;
        FigHeight = Ps.FigMaxHeight;
    end
    if AxesHeight < LegendTextHeight
        LegendScaleFactor = (AxesHeight - length(CDR3legend)/72) / TotalLegendHeight;
        Ps.LegendFontSize = floor(Ps.LegendFontSize * LegendScaleFactor);
        if Ps.LegendFontSize < 4 %Otherwise it's too small
            Ps.LegendFontSize = 4;
        end
    end
    
    resizeFigure(Gx, 'FigWidth', -1, 'FigHeight', FigHeight, 'Recenter', 'n');
    resizeSubplots(Gx, 'FigSpacer', 0.01);

    %----------------------------------------------------------------------
    %Rescaling horizontal components depending on legend

    if ~isempty(UnqCDR3seq{1}) && strcmpi(Ps.Legend, 'y')
        %Figure out what the text width would be
        HorzSpacer = 0.01;
        TempText = text(1, 1, UnqCDR3seq{1}, 'FontName', 'Courier', 'FontSize', Ps.LegendFontSize, 'Units', 'inches'); %Use courier for even spacing. Include spacing at end.
        TextExt = get(TempText, 'Extent');
        TextWidth = TextExt(3) + HorzSpacer; %inches, adding 0.01 in right side spacer
        delete(TempText)

        %Need to adjust Xlim to make sure it doesn't cross into legend
        PlotPosition = get(Ax, 'position');
        Xlim = get(Ax, 'Xlim');
        if isempty(Ps.Xmax) || Ps.Xmax == 0
            NewXlim2 = Xlim(1) + (MaxX - Xlim(1)) * PlotPosition(3) / (PlotPosition(3) - TextWidth);
            if NewXlim2 > Xlim(2)
                Xlim(2) = ceil(NewXlim2/5)*5;
            end
        else %Override Xlim
            if Ps.Xmax < Xlim(1)
                Xlim(2) = Xlim(1) + 1;
            else
                Xlim(2) = Ps.Xmax;
            end
        end
           
        if Xlim(2) < 10
            Xlim(2) = 10;
        end
        if mod(Xlim(2), 5) == 0
            Xlim(2) = Xlim(2) + 2; %We do this to prevent xlabel from going out of bounds, which causes a plot size reshift.
        end
        Xlim(1) = 0 - 0.05*Xlim(2); %Ensures that 0 start consistent distance from left edge.

        %You want at most 10 values on X, and multiple of 5's if possible)
        Xincr = ceil(ceil(Xlim(2)/10)/5)*5;
        XTickVal = 0:Xincr:Xlim(2);
        XTickLab = cell(size(XTickVal));
        for w = 1:length(XTickLab)
            XTickLab{w} = num2str(XTickVal(w), '%d');
        end
        set(Ax, 'XLim', Xlim, ...
                'XTick', XTickVal, ...
                'XTickLabel', XTickLab, ...
                'XMinorTick', 'off');

        %XY coordinate, right-aligned. Use inch, with resect to plot botleft corner.
        PlotPosition = get(Ax, 'Position');
        XYcoor = zeros(size(UnqCDR3seq, 1), 2);
        XYcoor(:, 2) = PlotPosition(4) - LegendTextHeight*(1:length(UnqCDR3seq)); %Vertical anchor points
        XYcoor(:, 1) = PlotPosition(3) - HorzSpacer; %Horizontal anchor points, based from the right border of plot. With 0.01 in spacer.

        %Add the legend text on the plots
        for j = 1:length(UnqCDR3seq)
            TextName = strrep(CDR3legend{j}, '_', '\_'); %Ensures any underscores prevent making subscripts
            text(XYcoor(j, 1), XYcoor(j, 2), TextName, 'FontWeight', 'bold', 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'bottom', 'FontName', 'Courier', 'FontSize', Ps.LegendFontSize, 'Color', UnqDotColor(j, :), 'Units', 'inch');
        end
    end
    
    SavePre = sprintf('.Grp%d', UnqGrpNum(y));
    FullSaveName = prepSaveTarget('SaveAs', Ps.SaveAs, 'SavePrefix', SavePre);
    FullSaveName = savePlot(Gx, 'SaveAs', FullSaveName); %300dpi, png.

    ImageGrpNums(y) = UnqGrpNum(y);
    ImageNames{y} = FullSaveName;
end
%Close plot if set to invisible
if strcmpi(Ps.Visible, 'off') 
    close(Gx);
end

%Create the tree search table
DelLoc = ImageGrpNums == 0;
ImageGrpNums(DelLoc) = [];
ImageNames(DelLoc) = [];

TreeSearchTable = getTreeSearchTable(VDJdata, VDJheader);
[~, ~, TableIdx] = intersect(ImageGrpNums, cell2mat(TreeSearchTable(2:end, 1)));
TreeSearchTable = TreeSearchTable([1; TableIdx+1], :);

%Remove the NewFilePath from the TreeFiles
for j = 1:length(ImageNames)
    [FilePath, FileName, ~] = parseFileName(ImageNames{j});
    ImageNames{j} = FileName;
end
TreeSearchTable(2:end, end) = ImageNames;
writeDlmFile(TreeSearchTable, [FilePath 'TreeSearchTable.csv']);

varargout{1} = ImageGrpNums;
varargout{2} = ImageNames;
varargout{3} = TreeSearchTable;
