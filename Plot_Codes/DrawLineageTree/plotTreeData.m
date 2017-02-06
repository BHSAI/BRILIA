%plotTreeData will plot the lineage trees. Requires that the TreeData is
%obtained from VDJdata, using getTreeData function. To use a custom color
%scheme for each tree node, assemble the Mx3 RGB matrix separately and then
%use that as a DotColor input.
%
%  [Gxs, Axs] = plotTreeData
%
%  [Gxs, Axs] = plotTreeData(FileName)
%
%  [Gxs, Axs] = plotTreeData(TreeData,TreeHeader)
%
%  [Gxs, Axs] = plotTreeData(TreeData,TreeHeader,Param,Value,...)
%
%  [Gxs, Axs] = plotTreeData(Param,Value,...)
%
%  [Gxs, Axs] = plotTreeData(...,P)
% 
%  P = plotTreeData('getinput)
%
%  INPUT
%    TreeData: cell containing lineage info. If not provided, will ask user
%      to select the VDJdata file returned by BRILIA (should be a CSV file)
%    TreeHeader: cell containing data header information
%    P: a structure containing P.(Param) = Value information. This input
%       option might be easier for adjusting plot parameters. You can also
%       use the direct param, cell pair input listed below. To make a
%       default P structure, use P = plotTreeData('GetInput');
%
%    Param-Value pair can be inputs, and below explains what they do.
%
%    *You can filter what trees to plot using these param-value pairs: 
%      Parmaeter Name  Value       Description
%      --------------- ---------   ----------------------------------------
%      GetGrpNum       [GrpNums]   Groups to plot. If empty, default is top
%                                    10 largest groups.
%      GetSeqNum       [SeqNums]   Groups with these SeqNums to plot. If
%                                    empty, default is top 10 largest
%                                    groups.
%      GetSizeRange    [Min,Max]   Groups with min to max # of seq to plot.
%                                    Default is range that gives you top 10
%                                    largest sizes.
%      GetCDR3seq      {'CDR3seq'} Groups with this CDR3 seq to plot. Can
%                                    be a char or cell of strings. Default
%                                    is {''}, do not filter.
%
%    *You can adjust the tree drawing using these param-value pairs:
%      Parmaeter Name  Value       Description
%      --------------- ---------   ----------------------------------------
%      DistanceUnit    'shm','ham' Show node-to-node distance either in SHM
%                                    or hamming distance
%      DotMaxSize      N >= 100    Max dot area size, in (pixel^2). Default
%                                    is 800.
%      DotScalor       N >= 1      Multiply template count by this to get
%                                    DotSize per each sequence. Will not
%                                    exceed DotMax. Default is 30.
%      DotColorMap     Mx3 matrix  A RGB colormap matrix used to decide how
%                                    to color each tree node. Default is a 
%                                    jet colormap.
%      Legend          'y','n'     Will draw a color-coded CDR3 legend.
%                                    Default is y.
%      LendFontSize    10          Font size of the Legend
%      Sort            'y','n'     Will sort the tree to cluster related
%                                    branches together.
%      TreeStyle       'triangle'  Draw trees with triangle edges. Default.
%                      'square'    Draw trees with square edges
%      Xmax            N           SHM or HAM distance limit for the X
%                                    axis. Default is 20.
%      Yincr           N           Vertical spacing between nodes. Default
%                                    is 0.125.
%
%    *Other general plotting param-value pairs:
%      Parmaeter Name  Value       Description
%      --------------- ---------   ----------------------------------------
%      FontSize        10          Font size of X and Y axes
%      FontName        Arial       Font name of X and Y axes
%      FigWidth        inches      Width of the whole figure
%      FigMaxHeight    inches      Maximum height of the whole figure
%
%    *If you want to save, use these param-value paris:
%      Parmaeter Name  Value       Description
%      --------------- ---------   ----------------------------------------
%      Save            'n','y'     To save or not to save. Default 'n'.
%      SaveName        String      The file name to save everything. The
%                                    file name will append the GrpNum to
%                                    it. EX: if SaveName is 'A', then the
%                                    files will be saved as A.Grp1.tif. If
%                                    empty and Save = 'y', will ask user to
%                                    select folder and file name.
%      Format          'tif','png' Image format type to save the figure as.
%                      'jpg','fig'   Default is 'tif'.
%      DPI             N           Dots per square inch. Default is 300.
%
%  OUTPUT
%    Gxs: figure handle(s) in a cell matrix
%    Axs: axes handle(s) in a cell matrix
%
%  NOTE
%    The filter param-value pairs uses an OR logic with each other,
%    selecting any group that matches any of the filter rules. Example:
%    GrpNum = 10 OR SeqNum = 12 will plot 1 OR 2 trees, depending on if
%    Seq 12 is or isn't within Grp 10.
%
%  EXAMPLE 
%    Go to the folder BRILIA/Examples_Files/BRILIA/ and try the following:
%      FileName = 'Ex4_SimMouseBCR_FullLength.BRILIAv2.0.4.csv';
%      P = plotTreeData('getinput');   %Get the default input struct only
%      P.GetGrpNum = 1;                %Plot tree for Group #1
%      P.TreeStyle = 'square';         %Draw tree with 90 degree lines
%      P.DotColorMap = jet/2;          %Use dark, jet colormap scheme
%      [Gxs,Axs] = plotTreeData(FileName,P); %Plots the tree using input P
%
%    To save the plots as they are made, modify P as such:
%      P.Save = 'y';                   %Yes, save the file.
%      P.Format = 'tif';               %Select file format
%      P.DPI = 400;                    %Select resolution
%      P.SaveName = 'Sample';          %A file name prefix
%      [Gxs,Axs] = plotTreeData(FileName,P);
%
%  See also getTreeData, exportFigure
function varargout = plotTreeData(varargin)
%The special case check for extracting input to plotTreeData
JustGettingInput = 0;
if nargin == 1 && ischar(varargin{1})
    if strcmpi(varargin{1},'getinput')
        JustGettingInput = 1;
        varargin = {}; %Want to get defaults.
    end
end

if JustGettingInput == 0
    %See if user gave TreeData
    if isempty(varargin) %Need to find file
        [TreeData,TreeHeader] = getTreeData;
    else
        if isempty(varargin{1}) %Empty file name, or TreeData and/or TreeHeader was given.
            [TreeData,TreeHeader] = getTreeData;
            varargin(1) = [];
            if mod(length(varargin),2) ~= 0 && ~isstruct(varargin{1}) %Maybe user gave empy TreeHeader too
                varargin(1) = []; %Delete to make varargin even again
            end
        elseif isstruct(varargin{1}) %A structre input without file name
            [TreeData,TreeHeader] = getTreeData;            
        elseif ischar(varargin{1}) && ~isempty(regexp(varargin{1},'\.','once')) %Filename was given with a period extension
            [TreeData,TreeHeader] = getTreeData(varargin{1});
            varargin(1) = [];
        elseif ischar(varargin{1}) &&  isempty(regexp(varargin{1},'\.','once')) %Parameter name was given first instead without TreeData or FileName
            [TreeData,TreeHeader] = getTreeData;
        elseif nargin >= 2 && iscell(varargin{1}) && iscell(varargin{2}) %TreeData and TreeHeader were given
            TreeData = varargin{1};
            TreeHeader = varargin{2};
            varargin(1:2) = [];
        else
            error('getTreeData: Check the inputs');
        end
    end
    TH = getTreeHeaderVar(TreeHeader); %Warning! Do not use getTreeHeaderVar AND H = getHeaderVar(VDJheader);

    %If structure param-value pairs are givert, convert to cell
    if ~isempty(varargin) && isstruct(varargin{1})
        S = varargin{1};
        Svalue = struct2cell(S)';
        Sparam = fieldnames(S)';
        Scell = [Sparam; Svalue];
        varargin = Scell(:)';
    end
end

%Parse the input
P = inputParser;
%Filter parameters (to draw a subset of trees)
addParameter(P,'GetGrpNum',[],@(x) isempty(x) || isnumeric(x));
addParameter(P,'GetSeqNum',[],@(x) isempty(x) || isnumeric(x));
addParameter(P,'GetSizeRange',[],@(x) isempty(x) || isnumeric(x));
addParameter(P,'GetCDR3seq',[],@(x) isempty(x) || ischar(x) || iscell(x));
%Display parameters (that directly affects tree)
addParameter(P,'DistanceUnit','shm',@(x) any(validatestring(lower(x),{'shm','ham'})));
addParameter(P,'DotMaxSize',800,@(x) isnumeric(x) && x >= 1);
addParameter(P,'DotScalor',30,@(x) isnumeric(x) && x >= 1);
addParameter(P,'DotColorMap',[],@(x) isempty(x) || (isnumeric(x) && size(x,2) == 3));
addParameter(P,'Legend','y',@(x) any(validatestring(lower(x),{'y','n'})));
addParameter(P,'LegendFontSize',10,@isnumeric)
addParameter(P,'Sort','y',@(x) any(validatestring(lower(x),{'y','n'})));
addParameter(P,'TreeStyle','triangle',@(x) any(validatestring(lower(x),{'triangle','square'})));
addParameter(P,'Xmax',0,@isnumeric);
addParameter(P,'Yincr',0.125,@isnumeric);
%General plotting parameters (that doesn't need tree)
addParameter(P,'FontSize',10,@isnumeric);
addParameter(P,'FontName','Arial',@ischar);
addParameter(P,'FigWidth',3.3,@isnumeric);
addParameter(P,'FigMaxHeight',5,@isnumeric);
%Saving parameters
addParameter(P,'Save','n',@(x) any(validatestring(lower(x),{'y','n'})));
addParameter(P,'SaveName','',@ischar);
addParameter(P,'Format','tif',@(x) any(validatestring(lower(x),{'tif','jpg','png','fig'})));
addParameter(P,'DPI',300,@isnumeric);

parse(P,varargin{:});
P = P.Results; %Remove Results field

if JustGettingInput == 1
    varargout{1} = P;
    return;
end

%==========================================================================
%Filter tree data, returning grp numbers to evaluate

SeqNum = cell2mat(TreeData(:,TH.SeqNumLoc));
GrpNum = cell2mat(TreeData(:,TH.GrpNumLoc));
GrpSize = zeros(size(GrpNum));
CDR3seq = TreeData(:,TH.CDR3seqLoc);

%Determine grp sizes per each entity
UnqGrpNum = unique(GrpNum);
for y = 1:length(UnqGrpNum)
    GrpIdx = UnqGrpNum(y) == GrpNum;
    GrpSize(GrpIdx) = sum(GrpIdx);
end

%Determine if you need to find top 10, or specified range
if isempty(P.GetSizeRange) && isempty(P.GetGrpNum) && isempty(P.GetSeqNum) %Find Top 10
    UnqGrpSize = sort(unique(GrpSize),'descend'); %Max size is 1st
    if length(UnqGrpSize) > 10
        MinGrpSize = UnqGrpSize(10);
    else
        MinGrpSize = UnqGrpSize(end);
    end
    ValidGrpSize = GrpSize >= MinGrpSize;
elseif ~isempty(P.GetSizeRange) %Use specified size range
    ValidGrpSize = GrpSize >= P.GetSizeRange(1) & GrpSize <= P.GetSizeRange(end);
else
    ValidGrpSize = zeros(size(TreeData,1),1,'logical');
end

%Determine what groups are valid, based on grp num
if isempty(P.GetGrpNum)
    ValidGrpNum = ValidGrpSize;
else
    ValidGrpNum = zeros(size(GrpNum),'logical');
    for y = 1:length(P.GetGrpNum)
        ValidGrpNum = ValidGrpNum | GrpNum == P.GetGrpNum(y);
    end
end

%Determine what groups are valid, based on seq num
if isempty(P.GetSeqNum)
    ValidSeqNum = ValidGrpSize;
else
    ValidSeqNum = zeros(size(SeqNum),'logical');
    [~,GetThese,~] = intersect(SeqNum,P.GetSeqNum);
    
    %Find all group members that have these sequences
    GetGrpNum = GrpNum(GetThese);
    for y = 1:length(GetGrpNum)
        GrpIdx = GetGrpNum(y) == GrpNum;
        ValidSeqNum(GrpIdx) = 1;
    end
end

%Determine what groups are valid, based on CDR3seq
if isempty(P.GetCDR3seq)
    ValidCDR3seq = zeros(size(TreeData,1),1,'logical');
else
    ValidCDR3seq = zeros(size(CDR3seq),'logical');
    for y = 1:length(UnqGrpNum)
        GrpIdx = UnqGrpNum(y) == GrpNum;
        UnqGrpCDR3seq = unique(CDR3seq(GrpIdx));
        HaveCDR3seq = findCell(UnqGrpCDR3seq,P.GetCDR3seq);
        if max(HaveCDR3seq) > 0
            ValidCDR3seq(GrpIdx) = 1;
        end
    end
end
    
%Determine what groups numbers to keep
KeepThese = ValidSeqNum | ValidGrpNum | ValidGrpSize | ValidCDR3seq;
EvalGrpNum = unique(GrpNum(KeepThese));

if isempty(EvalGrpNum); 
    disp('No trees fit the filter criteria.');
    return
end

%End of filtering

%==========================================================================
%Process each tree, group by group
SaveNamePre = ''; %Prefix for the file name to save. Initialized with ''.
Axs = cell(length(EvalGrpNum),1);
Gxs = cell(length(EvalGrpNum),1);
for y = 1:length(EvalGrpNum)
    %----------------------------------------------------------------------
    %Determining plotting info, coord, size, etc

    %Isolate group data and name
    GrpIdx = GrpNum == EvalGrpNum(y);
    TreeDataT = TreeData(GrpIdx,:);
    GrpName = TreeDataT{1,TH.GrpNameLoc};

    %Obtain the X-Y coordinates for plotting the tree
    if strcmpi(P.DistanceUnit,'shm')
        DistLoc = TH.SHMdistLoc;
    else %defaults to hamming distance
        DistLoc = TH.HAMdistLoc;
    end
    AncMap = cell2mat(TreeDataT(:,[TH.ChildSeqNumLoc TH.ParentSeqNumLoc DistLoc])); %Simple ancestry map matrix
    AncMap = renumberAncMap(AncMap); %Ancestry maps to relative position now
    if lower(P.Sort(1)) == 'y' %Sort helps to group related sequences
        [AncMap,SortIdx] = sortrows(AncMap,[2 3 1]);
        TreeDataT = TreeDataT(SortIdx,:);
        AncMap = renumberAncMap(AncMap);
    end
    TreeCoord = calcTreeCoord(AncMap);
    
    %Make the legend text and get default color scheme using CDR3
    CDR3seqT = TreeDataT(:,TH.CDR3seqLoc);
    [CDR3legend, UnqCDR3seq] = makeTreeLegend_CDR3(CDR3seqT);
    
    %Map a unique CDR3 colors to TreeDataT based on the order UnqCDR3seq.
    if isempty(P.DotColorMap) %Make the default colors
        [DotColor, UnqDotColor] = mapDotColor_CDR3(CDR3seqT,UnqCDR3seq);
    else
        [DotColor, UnqDotColor] = mapDotColor_CDR3(CDR3seqT,UnqCDR3seq,'ColorMap',P.DotColorMap);
    end
    
    %Calc the dot sizes, resizing dots if one exceeds the max dot size.
    DotSize = cell2mat(TreeDataT(:,TH.TemplateLoc))*P.DotScalor;
    CurMaxDotSize = max(DotSize);
    if CurMaxDotSize > P.DotMaxSize %Need rescaling
        DotSize = DotSize * P.DotMaxSize / CurMaxDotSize;
    end    
    
    %----------------------------------------------------------------------
    %Drawing the lineage tree

    %Setup template figure
    [Gx,Ax] = formatPlotSize([],P.FigWidth,P.FigMaxHeight); %Creates axes with inches unit!
    Gxs{y} = Gx; %Save to output
    Axs{y} = Ax; %Save to output
    set(Ax,'XTickLabelMode','auto',...
           'YTickLabelMode','manual',...
           'TickLength',[0.005 0.005],...
           'TickDir','both',...
           'YTickLabel','',...
           'YTick',[],...
           'FontName',P.FontName,...
           'FontSize',P.FontSize,...
           'box','on');

    %Add the X axis and title labels
    if strcmpi(P.DistanceUnit,'shm')
        xlabel(Ax,'SHM Distance','FontName',P.FontName,'FontSize',P.FontSize);
    else
        xlabel(Ax,'HAM Distance','FontName',P.FontName,'FontSize',P.FontSize);
    end
    if ~isempty(GrpName)
        GrpName = strrep(GrpName,'_','\_'); %Prevent underscores from being subscript
        title(Ax,GrpName,'FontName',P.FontName,'FontSize',P.FontSize+1);
    end

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
        if strcmpi(P.TreeStyle,'square') %Uses bracket-like drawing
            plot([X0 X0 X1],[Y0 Y1 Y1],'k','LineWidth',1)
        elseif strcmpi(P.TreeStyle,'triangle') %Uses direct lines between nodes
            plot([X0 X1],[Y0 Y1],'k','LineWidth',1)
        end
    end

    %Draw tree nodes and leaves
    [~,IdxT] = sort(DotSize,'descend'); %Sort to prevent covering up smaller dot
    scatter(Ax,TreeCoord(IdxT,1),TreeCoord(IdxT,2),DotSize(IdxT,:),DotColor(IdxT,:),'fill')
    hold(Ax,'off')

    %Reset the labels and axis lengths
    MaxX = max(TreeCoord(:,1));
    MaxY = max(TreeCoord(:,2));
    set(Ax,'Ylim',[0 MaxY+1]);
    set(Ax,'Xlim',[-2,MaxX+2]);

    %----------------------------------------------------------------------
    %Preparing tree plot height to reduce wasted space
    VertSpacer = 1/72; %Vertical spacer between legend entries
    HorzSpacer = 0.04; %Horizontal spacer from right of text to right plot border
    
    %Determine height of legend, adding 1pt spacer and 1 extra full text
    TextHeight = P.LegendFontSize / 72 + VertSpacer; %Add 1pt spacer, convert pt to inch using 72pt = 1 in
    LegendHeight = (length(UnqCDR3seq)+1) * TextHeight; %Height of all legend
   
    %Determine height of tree
    if P.Yincr == 0 %Autoscale by font size
        P.Yincr = TextHeight;
    end
    TreeHeight = (MaxY+1) * P.Yincr; %Height of all tree
    
    %Determine required plot size, ensuring it does not exceed P.FigMaxHeight.
    ReqPlotHeight = max([TreeHeight LegendHeight]); %Required plot area height
    PlotBorder = get(Ax,'TightInset'); %Left Bot Right Top border in inches
    TotFigHeight = sum(PlotBorder([2 4])) + ReqPlotHeight; %See if the ReqPlotHeight makes figure taller than P.FigMaxHeight
    if TotFigHeight > P.FigMaxHeight %Need size rescaling
        TotFigHeight = P.FigMaxHeight;
        ReqPlotHeight = P.FigMaxHeight - sum(PlotBorder([2 4]));
        
        %If the plot height was adjusted due to legend, resize legend font size
        if LegendHeight > TreeHeight
            TextHeight = ReqPlotHeight / (length(UnqCDR3seq) +1);
            P.LegendFontSize = floor(TextHeight * 72 - VertSpacer);
            if P.LegendFontSize <= 4 %Don't make it smaller than 4, sinc you can't read it
                P.LegendFontSize = 4;
            end
        end
    end
    
    %Update the plot size
    PlotPosition = get(Ax,'Position'); %Axes position in inches
    PlotPosition(4) = ReqPlotHeight;
    set(Ax,'Position',PlotPosition);
    formatPlotSize(Ax,P.FigWidth,TotFigHeight); %Sets the final figure height and width 

    %----------------------------------------------------------------------
    %Rescale X limits such that dots do not overlap with the legend text
    
    if ~isempty(UnqCDR3seq{1}) && strcmpi(P.Legend,'y')
        %Figure out what the text width would be
        TempText = text(1,1,UnqCDR3seq{1},'FontName','Courier','FontSize',P.LegendFontSize,'Units','inches'); %Use courier for even spacing. Include spacing at end.
        TextExt = get(TempText,'Extent');
        TextWidth = TextExt(3) + HorzSpacer; %inches, adding 0.01 in right side spacer
        delete(TempText)

        %Need to adjust Xlim to make sure it doesn't cross into legend
        Xlim = get(Ax,'Xlim');
        if isempty(P.Xmax) || P.Xmax == 0
            NewXlim2 = Xlim(1) + (MaxX - Xlim(1)) * PlotPosition(3) / (PlotPosition(3) - TextWidth);
            if NewXlim2 > Xlim(2)
                Xlim(2) = ceil(NewXlim2/5)*5;
            end
        else %Override Xlim
            if P.Xmax < Xlim(1)
                Xlim(2) = Xlim(1) + 1;
            else
                Xlim(2) = P.Xmax;
            end
        end
        if mod(Xlim(2),5) == 0
            Xlim(2) = Xlim(2) + 2; %We do this to prevent xlabel from going out of bounds, which causes a plot size reshift.
        end

        %You want at most 10 values on X, and multiple of 5's if possible)
        Xincr = ceil(ceil(Xlim(2)/10)/5)*5;
        XTickVal = [0:Xincr:Xlim(2)];
        XTickLab = cell(size(XTickVal));
        for w = 1:length(XTickLab)
            XTickLab{w} = num2str(XTickVal(w),'%d');
        end
        set(Ax,'XLim',Xlim,...
               'XTick',XTickVal,...
               'XTickLabel',XTickLab,...
               'XMinorTick','off');
        formatPlotSize(Ax,P.FigWidth,TotFigHeight); %Ensure labels in X doesn't get cut off

        %XY coordinate, right-aligned. Use inch, with resect to plot botleft corner.
        PlotPosition = get(Ax,'Position');
        XYcoor = zeros(size(UnqCDR3seq,1),2);
        XYcoor(:,2) = PlotPosition(4) - TextHeight*[1:length(UnqCDR3seq)]; %Vertical anchor points
        XYcoor(:,1) = PlotPosition(3) - HorzSpacer; %Horizontal anchor points, based from the right border of plot. With 0.01 in spacer.

        %Add the legend text on the plots
        for j = 1:length(UnqCDR3seq)
            TextName = strrep(CDR3legend{j},'_','\_'); %Ensures any underscores prevent making subscripts
            text(XYcoor(j,1),XYcoor(j,2),TextName,'FontWeight','bold','HorizontalAlignment','Right','VerticalAlignment','middle','FontName','Courier','FontSize',P.LegendFontSize,'Color',UnqDotColor(j,:),'Units','inch');
        end
    end
    
    %Save the figure, if the user wants
    if strcmpi(P.Save,'y')
        if isempty(SaveNamePre) %Need to establish the prefix, save path, file ext here.
            if isempty(P.SaveName) %Need to ask users to select file name
                [SaveFile, SavePath] = uiputfile('*.tif;*.png;*.jpg;*.fig');
                [SavePath, SaveFile, SaveExt] = parseFileName([SavePath SaveFile]);
            else
                [SavePath, SaveFile, SaveExt] = parseFileName(P.SaveName);
                if isempty(SaveExt)
                    SaveExt = ['.' P.Format]; %Use the format option specified in the input. P.Format cannot be empty.
                    if SaveExt(2) == '.'; SaveExt(2) = []; end %Prevents ..jpg file format nuissances
                end
            end
            DotLoc = find(SaveFile == '.');
            if isempty(DotLoc)
                SaveNamePre = SaveFile;
            else
                SaveNamePre = SaveFile(1:DotLoc(end)-1);
            end
        end
                
        %Assemble the full save name, and save depending on file ext
        FullSaveName = sprintf('%s%s.Grp%d%s',SavePath,SaveNamePre,EvalGrpNum(y),SaveExt);
        switch SaveExt
            case '.tif'
                print(Gx,FullSaveName,'-dtiff',['-r' num2str(P.DPI)]);
            case '.jpg'
                print(Gx,FullSaveName,'-djpeg',['-r' num2str(P.DPI)]);
            case '.png'
                print(Gx,FullSaveName,'-dpng',['-r' num2str(P.DPI)]);                
            case '.fig'
                saveas(Gx,FullSaveName);
        end
    end    
end

if nargout >= 1
    varargout{1} = Gxs; %Figure handles
    if nargout >= 2
        varargout{2} = Axs; %Axes handles
    end
end
