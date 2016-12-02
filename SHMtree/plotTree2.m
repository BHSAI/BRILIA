%plotTree will take in a ancestral tree matrix, and then plot the tree.

%  [Gx,Ax] = plotTree2(AncMap,TreeName,CDR3seq,SortMode,CDR3clrCell), where
%  AncMap = the Nx3 matrix created by buildTreeLink. TreeName is the text
%  of the tree name (recommend to put something descriptive), CDR3seq is a
%  cell array of text for the CDR3, and SortMode is 'sort','none' if you
%  want to sort the tree based on similar CDR3.

function [Gx,Ax] = plotTree2(AncMap,varargin)
%==========================================================================
%Parse the input
P = inputParser;
addRequired(P,'AncMap',@isnumeric);
addOptional(P,'TreeName','',@ischar);
addOptional(P,'CDR3seq','',@iscell);
addOptional(P,'DotClr',[],@isnumeric);
addParameter(P,'SortMode','sort',@(x) any(validatestring(x,{'sort','none'})));
addParameter(P,'FigWidth',3.3,@isnumeric);
addParameter(P,'FigMaxHeight',8,@isnumeric);
addParameter(P,'FigFontSize',9,@isnumeric);
addParameter(P,'SetXlim',0,@isnumeric); 
addParameter(P,'SetYincr',0,@isnumeric); %This will set the spacing between y's
addParameter(P,'ScaleFactor',50,@isnumeric); %This scales the dot size x ScaleFactor pixel^2.
addParameter(P,'MaxDotArea',800,@isnumeric); %This sets the maximum dot area, rescaling all dots. This overides ScaleFactor.
parse(P,AncMap,varargin{:});

AncMap = P.Results.AncMap;
TreeName = P.Results.TreeName;
CDR3seq = P.Results.CDR3seq;
DotClr = P.Results.DotClr;
SortMode = P.Results.SortMode;
FigWidth = P.Results.FigWidth;
FigMaxHeight = P.Results.FigMaxHeight;
FigFontSize = P.Results.FigFontSize;
SetXlim = P.Results.SetXlim;
SetYincr = P.Results.SetYincr;
ScaleFactor = P.Results.ScaleFactor;
MaxDotArea = P.Results.MaxDotArea;

%==========================================================================
%Getting information to build the tree

%Make sure AncMap has template count info. If not, add 1's to col 4.
if size(AncMap,2) <= 3 %If template count is provided
    AncMap = [AncMap ones(size(AncMap,1),1)];
end

%Sort and label unique CDR3seq. Assign dot color to for each CDR3seq.
[UnqCDR3,~,UnqCDR3align] = plotTree_getUnqCDR3(CDR3seq,AncMap);
if isempty(DotClr) %Create one
    [DotClr,RefClr] = plotTree_getDotClr(CDR3seq,UnqCDR3);
else %Dot color is already specified. Need to create RefClr for legend.
    RefClr = zeros(size(UnqCDR3,1),3);
    for k = 1:size(UnqCDR3,1)
        for j = 1:size(CDR3seq,1)
            if strcmpi(CDR3seq{j},UnqCDR3{k})
                RefClr(k,:) = DotClr(j,:);
                break
            end
        end
    end
end

%Set the dot area scaling factor so that dots aren't too big/small.
CurMaxDot = max(AncMap(:,4))*ScaleFactor;
if CurMaxDot > MaxDotArea %Need to rescale the ScaleFactor;
    ScaleFactor = ScaleFactor*MaxDotArea/CurMaxDot;
end

%Renumber AncMap for convenience, and to ensure root seq is near top
if strcmpi(SortMode,'sort')
    [AncMap,SortIdx] = sortrows(AncMap,[2 3 1]);
    DotClr = DotClr(SortIdx,:);
    AncMap = renumberAncMap(AncMap);
end
TreeCoord = calcTreeCoord(AncMap);

%==========================================================================
%Drawing the lineage tree

%Setup template figure
[Gx,Ax] = formatPlotSize([],FigWidth,FigMaxHeight); %Creates axes with inches unit!
set(Ax,'XTickLabelMode','auto',...
       'YTickLabelMode','manual',...
       'TickLength',[0.005 0.005],...
       'TickDir','both',...
       'YTickLabel','',...
       'YTick',[],...
       'FontName','Arial',...
       'FontSize',FigFontSize,...
       'box','on');

%Add the X axis and title labels
if ~isempty(TreeName)
    TreeName = strrep(TreeName,'_','\_'); %Prevent underscores from being subscript
    title(Ax,TreeName,'FontName','Arial','FontSize',FigFontSize+1);
end
xlabel(Ax,'Lineage Distance','FontName','Arial','FontSize',FigFontSize);

%Draw the lineage lines only
hold(Ax,'on')
for j = 1:size(AncMap,1)
    ParentLoc = AncMap(j,2);
    if ParentLoc == 0 %if it's the root, use TreeCoord value directly.
        ParentLoc = j;
    end
    X0 = TreeCoord(ParentLoc,1);
    Y0 = TreeCoord(ParentLoc,2);
    X1 = TreeCoord(j,1);
    Y1 = TreeCoord(j,2);
    plot([X0 X0 X1],[Y0 Y1 Y1],'k','LineWidth',1)
end

%Sort the variables by dot size to prevent cover ups. Draw tree nodes.
[~,SortIdxT] = sort(AncMap(:,4),'descend');
TreeCoordT = TreeCoord(SortIdxT,:);
AncMapT = AncMap(SortIdxT,:);
DotClrT = DotClr(SortIdxT,:);
scatter(Ax,TreeCoordT(:,1),TreeCoordT(:,2),AncMapT(:,4)*ScaleFactor,DotClrT,'fill')
hold(Ax,'off')

%Reset the labels and axis lengths
MaxX = max(TreeCoord(:,1));
MaxY = max(TreeCoord(:,2));
set(Ax,'Ylim',[0 MaxY+1]);
set(Ax,'Xlim',[-2,MaxX+1]);

%==========================================================================
%Preparing tree plot height to reduce wasted space

%Determine spacing between tree Y axis, in inches
TextHeight = (FigFontSize+1) * 0.014; %Point to inch conversion
if SetYincr == 0 %Autoscale by font size
    SetYincr = TextHeight;
end

%Determine the required plot axis height
TreeHeight = (MaxY+1) * SetYincr;
LegendHeight = (size(UnqCDR3,1)+1) * TextHeight;
ReqHeight = max([TreeHeight LegendHeight]);

%Determine required plot size, ensuring it does not exceed FigMaxHeight.
PlotPosition = get(Ax,'Position'); %Axes position in inches
PlotBorder = get(Ax,'TightInset'); %Left Bot Right Top border in inches
TotFigHeight = sum(PlotBorder([2 4])) + ReqHeight; %See if the reqheight makes figure taller than 10 inch, max
if TotFigHeight > FigMaxHeight
    ReqHeight = FigMaxHeight - sum(PlotBorder([2 4]));
    TotFigHeight = FigMaxHeight;
end
PlotPosition(4) = ReqHeight;
set(Ax,'Position',PlotPosition);
formatPlotSize(Ax,FigWidth,TotFigHeight);

%==========================================================================
%Rescale X such that it does not go into legend

if ~isempty(UnqCDR3{1})
    %Figure out what the text width would be
    TempStr = [CDR3seq{1} ' ']; %want a space at the end
    TempText = text(1,1,TempStr,'FontName','Courier','FontSize',FigFontSize); %Use courier for even spacing
    set(TempText,'Units','inches')
    TextExtInch = get(TempText,'Extent');
    delete(TempText)

    %Figure out the xlimits
    Xlimits = get(Ax,'Xlim');
    if SetXlim == 0   %Need to rescale Xlim to make sure it doesn't cross into legend        
        InchOvrData = (PlotPosition(3) - TextExtInch(3)) / diff(Xlimits); %You want this inch/data scaling
        Xlimits(2) = ceil(PlotPosition(3)/InchOvrData - MaxX); %use MaxX since you want the lax point to not exceed the beginning of text legene
        Xlimits(2) = ceil(Xlimits(2)/5)*5;    
        if Xlimits(2) < 20 %Just to standardize plot X sizes somewhat
            Xlimits(2) = 20;
        end
    elseif SetXlim > 0   %Override Xlimits
        if SetXlim ~= 0
            Xlimits(2) = SetXlim;
        end
    end

    %You want at most 10 values on X, and multiple of 5's if possible)
    Xincr = ceil(ceil(Xlimits(2)/10)/5)*5;
    XTickVal = [0:Xincr:Xlimits(2)];
    XTickLab = cell(size(XTickVal));
    for w = 1:length(XTickLab)
        XTickLab{w} = num2str(XTickVal(w),'%d');
    end
    set(Ax,'XLim',Xlimits,...
           'XTick',XTickVal,...
           'XTickLabel',XTickLab,...
           'XMinorTick','off');
    formatPlotSize(Ax,FigWidth,TotFigHeight); %Ensure labels in X doesn't get cut off
    
    %XY coordinate, right-aligned. Use inch, with resect to plot botleft corner.
    XYcoor = zeros(size(UnqCDR3,1),2);
    PlotSize = get(Ax,'Position');
    PlotHeight = PlotSize(4);
    Yincr = TextHeight;
    ReqLegHeight = Yincr * (size(UnqCDR3align,1)+1);
    if ReqLegHeight > PlotHeight
        Yincr = PlotHeight/(size(UnqCDR3align,1)+1);
    end
    XYcoor(:,2) = PlotSize(4) - Yincr*[1:length(UnqCDR3)];
    XYcoor(:,1) = PlotSize(3);
        
    %Add the text on the plots
    for j = 1:length(UnqCDR3)
        TextName = sprintf('%s ',UnqCDR3align(j,:));
        TextName = strrep(TextName,'_','\_');        
        text(XYcoor(j,1),XYcoor(j,2),TextName,'FontWeight','bold','HorizontalAlignment','Right','VerticalAlignment','middle','FontName','Courier','FontSize',FigFontSize,'Color',RefClr(j,:),'Units','inch');
    end
end