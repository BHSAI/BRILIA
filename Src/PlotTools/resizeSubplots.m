%resizeSubplots will resize all subplots such that each plot has the same
%and largest axes size, without running axes labels off the figure.
%
%  resizeSubplots(Gx, varargin)
%
%  INPUT
%    Gx: figure handle to modify
%    varargin: cell array of Param-Value pairs.
%    
%    Param           Value         Description
%    ----------      ---------     ----------------------------------------
%    ScaleVertical   'y' 'n'       If 'y', will scale subplots vertically
%                                    based on original plot outer heights. 
%                                  If 'n', will divide plot heights evenly.
%    ScaleHorizontal 'y' 'n'       If 'y', will scale subplot horizontally
%                                    based on original plot outer widths. 
%                                  If 'n', will divide plot width evenly.
%    MatchWidth      'y' 'n'       If 'y', will set the width of each plot
%                                    the same per each column.
%                                  If 'n', will just fill the column width.
%    HorzSpacer      Number [0, 1] Will add a horizontal spacer between
%                                    subplot plots
%    VertSpacer      Number [0, 1] Will add a veritcal spacer between
%                                    subplot plots
%    FigSpacer       Number [0, 1] Will add a perimiter spacer around all
%                                    subplots 
%
%  EXAMPLE
%    Gx = figure();
%    subplot(2, 2, 1)
%    set(gca, 'position',  [0 0.5 0.3 0.2]);
%    plot(1, 1, 'r')
%    subplot(2, 2, 2)
%    plot(1, 2, 'b')
%    subplot(2, 2, 4)
%    plot(3, 3, 'k')
%    ylabel('blah')
%    resizeSubplots(Gx, 'ScaleVertical', 'n', 'MatchWidth', 'y', 'FigSpacer', 0.01)

function varargout = resizeSubplots(varargin)
[Gx, varargin] = getOpenGCF(varargin{:});
P = inputParser;
addParameter(P, 'ScaleVertical', 'y', @(x) validatestring(lower(x), {'y', 'n', 'yes', 'no'}));
addParameter(P, 'ScaleHorizontal', 'y', @(x) validatestring(lower(x), {'y', 'n', 'yes', 'no'}));
addParameter(P, 'MatchWidth', 'y', @(x) validatestring(lower(x), {'y', 'n', 'yes', 'no'}));
addParameter(P, 'HorzSpacer', 0, @(x) isnumeric(x));
addParameter(P, 'VertSpacer', 0, @(x) isnumeric(x));
addParameter(P, 'FigSpacer', 0, @(x) isnumeric(x));

[Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInputT(P, varargin{:});
if ReturnThis
   varargout = {Ps, Pu, ExpPs, ExpPu};
   return;
end
ScaleVertical = Ps.ScaleVertical(1);
MatchWidth = Ps.MatchWidth(1);
HorzSpacer = Ps.HorzSpacer;
VertSpacer = Ps.VertSpacer;
FigSpacer = Ps.FigSpacer;

%Save prior units, xlim, ylim as these could change upon scaling.
Axs = findobj(Gx, 'type', 'axes');
PriorUnits = cell(length(Axs), 1);
PriorXLims = cell(length(Axs), 1);
PriorYLims = cell(length(Axs), 1);
for j = 1:length(Axs)
    PriorXLims{j} = get(Axs(j), 'XLim');
    PriorYLims{j} = get(Axs(j), 'YLim');
    PriorUnits{j} = get(Axs(j), 'Units');
    set(Axs(j), 'Units', 'normalized'); %Do everything in normalized units.
end

%In case there are overlapped axes, operate per unique overlapped axes.
[UnqAxCell, AllPositions, AllTights] = getUniqueAxes(Gx);

%Determine subplot grid coord.
%Going top to bot is row 1, 2, ... 
%Going left to right is column 1, 2, ... 

%Find the columns first
GridCoord = zeros(size(AllPositions, 1), 2);
XPos = [AllPositions(:, 1) sum(AllPositions(:, [1 3]), 2)];
Active = ones(size(XPos, 1), 1, 'logical');
MinPos = min(XPos(:,2));
CurCol = 1;
while max(Active) == 1
    CurRowIdx = XPos(:,1) < MinPos & Active;
    GridCoord(CurRowIdx, 2) = CurCol;
    CurCol = CurCol + 1;
    Active(CurRowIdx) = 0;
    MinPos = min(XPos(Active, 2));
end

%Find the rows next
YPos = [AllPositions(:, 2) sum(AllPositions(:, [2 4]), 2)];
Active = ones(size(YPos, 1), 1, 'logical');
MaxPos = max(YPos(:, 1));
CurCol = 1;
while max(Active) == 1
    CurColIdx = YPos(:, 2) > MaxPos & Active;
    GridCoord(CurColIdx, 1) = CurCol;
    CurCol = CurCol + 1;
    Active(CurColIdx) = 0;
    MaxPos = max(YPos(Active, 1));
end

%Calculate the outer positions required of subplots and getMaxPlotPosition
AllOuterPos = zeros(size(AllPositions, 1), 4);
AllOuterPos(:, 3) = (1 - 2*FigSpacer) / max(GridCoord(:, 2)); %Determine W as equal widths

%Determine H of outer positions 
if strcmpi(ScaleVertical, 'y') %scaled by subplot heights
    for c = 1:max(GridCoord(:, 2))
        Idx = GridCoord(:, 2) == c;
        Vert = AllPositions(Idx, 4) + sum(AllTights(Idx, [2 4]), 2); %Vertical heights, including tight insets
        AllOuterPos(Idx, 4) = (1 - 2*FigSpacer) * (Vert / sum(Vert)); 
    end
else %scaled equally by number of subplots
    for c = 1:max(GridCoord(:, 2))
        Idx = GridCoord(:, 2) == c;
        AllOuterPos(Idx, 4) = (1 - 2*FigSpacer) / sum(Idx);
    end
end

%Determine the X and Y for all outer positions
for c = 1:max(GridCoord(:, 2))
    Idx = GridCoord(:, 2) == c; %Same column
    [~, SortIdx] = sort(GridCoord(Idx, 1)); %Get sort order by row number
    AllOuterPosT = AllOuterPos(Idx, :); %Get outer Pos for this column
    AllOuterPosT = AllOuterPosT(SortIdx, :);%Sort by row order
    AllOuterPosT(:, 1) = AllOuterPosT(:, 3) * (c - 1) + FigSpacer; %Determine X pos
    AllOuterPosT(:, 2) = (1 - FigSpacer) - cumsum(AllOuterPosT(:, 4)); %Determine Y pos
    [~, RevSortIdx] = sort(SortIdx); %Undo the sort by row
    AllOuterPos(Idx, :) = AllOuterPosT(RevSortIdx, :); %Replace with new values
end

%Rescale W to account for HorzSpacer
NumCol = max(GridCoord(:, 2));
RescaleW = ((1 - 2*FigSpacer) - HorzSpacer*(NumCol - 1)) / (1 - 2*FigSpacer);
AllOuterPos(:, 3) = AllOuterPos(:, 3)*RescaleW; %RescaleW for HorzSpacer
AdjIdx = GridCoord(:, 2) > 1;
AllOuterPos(AdjIdx, 1) = AllOuterPos(AdjIdx, 1) + HorzSpacer/2;

%Rescale H to account for VertSpacer
for c = 1:max(GridCoord(:, 2))
    Idx = find(GridCoord(:, 2) == c);
    NumRow = length(Idx);
    RescaleH = ((1 - 2*FigSpacer) - VertSpacer*(NumRow - 1)) / (1 - 2*FigSpacer);
    AllOuterPos(Idx, 4) = AllOuterPos(Idx, 4)*RescaleH; %RescaleH for VertSpacer
    AdjIdx = GridCoord(Idx, 1) ~= max(GridCoord(Idx, 1));
    AllOuterPos(Idx(AdjIdx), 2) = AllOuterPos(Idx(AdjIdx), 2) + VertSpacer/2;
end

%Determine the coordinates for each subplot that'll fill
AllMaxPos = zeros(size(AllOuterPos, 1), 4);
for j = 1:size(AllOuterPos, 1)
    AllMaxPos(j, 1) = AllOuterPos(j, 1) + AllTights(j, 1);
    AllMaxPos(j, 2) = AllOuterPos(j, 2) + AllTights(j, 2);
    AllMaxPos(j, 3) = AllOuterPos(j, 3) - AllTights(j, 1) - AllTights(j, 3);
    AllMaxPos(j, 4) = AllOuterPos(j, 4) - AllTights(j, 2) - AllTights(j, 4);
end

%Align horizontally within each column, if specified
if strcmpi(MatchWidth, 'y')
    for c = 1:max(GridCoord(:, 2))
        Idx = find(GridCoord(:, 2) == c);
        PosStart = max(AllMaxPos(Idx, 1)); %want the furthest start
        PosEnd = min(sum(AllMaxPos(Idx, [1 3]), 2)); %want the shorter end pos
        for j = 1:length(Idx)
            CurPosStart = AllMaxPos(Idx(j), 1);
            if CurPosStart ~= PosStart
                LeftTrim = PosStart - CurPosStart;
            else
                LeftTrim = 0;
            end

            CurPosEnd = sum(AllMaxPos(Idx(j), [1 3]), 2);
            if CurPosEnd ~= PosEnd
                RightTrim = CurPosEnd - PosEnd;
            else
                RightTrim = 0;
            end
            Xwidth = AllMaxPos(Idx(j), 3) - LeftTrim - RightTrim;

            AllMaxPos(Idx(j), 1) = PosStart;
            AllMaxPos(Idx(j), 3) = Xwidth;
        end
    end
end

%Resize all subplot
for j = 1:length(UnqAxCell)
    setAxes(UnqAxCell{j}, 'Position', AllMaxPos(j, :));
end

%Return the axes units, xlim, ylim to prior values
for j = 1:length(Axs)
    set(Axs(j), 'Units', PriorUnits{j}, 'XLim', PriorXLims{j}, 'YLim', PriorYLims{j});
end

