%labelSubplots will add the labels on the subplots at the top left corner
%of each plot. Can specify the labels to go left to right (horizontal) or
%top to bottom (vertical).
%
%  labelSubplots(Gx, param, value, ...)
%
%  labelSubplots(Gx, 'remove')
%
%  labelSubplots(Gx, 'shift', HShift, VShift)
%
%  INPUT
%    Gx: figure or axes handle
%    'remove': deletes all the panel labels
%    'shift': perform shifting only of the current panel labels
%    HShift: shifts labels horizontally. 
%    VShift: shifts labels vertically. 
%    
%    Param        Value (*)      Description
%    -----------  -------------  ------------------------------------------
%    Direction    * 'hor'        Labels a-z horizontally
%                   'ver'        Labels a-z vertically
%    FontName     * 'Arial'      Font name for labels
%    FontSize     * 10           Font size for labels
%    FontCase     * 'lower'      Lowercase panel labels
%                   'upper'      Uppercase panel labels
%    HShift       * 0            Shift labels left(-)/right(+)
%    VShift       * 0            Shift labels down(-)/up(+)
%
%  EXAMPLE
%    Gx = loadDemoPlot;
%    labelSubplots(Gx, 'FontName', 'Courier', 'Direction', 'ver');
%
function varargout = labelSubplots(varargin)
[Gx, varargin] = getOpenGCF(varargin{:});
if isempty(Gx); return; end

P = inputParser;
addParameter(P, 'Direction', 'hor',    @(x) startsWith(x, {'hor', 'ver'}, 'ignorecase', 1));
addParameter(P, 'FontName',  'Arial',  @ischar);
addParameter(P, 'FontSize',  10,       @(x) isnumeric(x) && x > 2);
addParameter(P, 'FontCase',  'lower',  @(x) ismember(lower(x), {'upper', 'lower'}));
addParameter(P, 'HShift',    0,        @isnumeric);
addParameter(P, 'VShift',    0,        @isnumeric);

%Check if this is for removing labels
Action = 0; %0 = add labels, 1 = delete labels, 2 = shift labels
for j = 1:length(varargin)
    if ischar(varargin{j}) && ismember(lower(varargin{j}), {'remove', 'delete', 'del', 'rem'})
        Action = 1;
        varargin = {}; 
        break
    elseif ischar(varargin{j}) && strcmpi(varargin{j}, 'shift')
        Action = 2;
        varargin = {};
        break
    end
end
[Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInput(P, varargin{:});
if ReturnThis
   varargout = {Ps, Pu, ExpPs, ExpPu};
   return
end
Direction = Ps.Direction;
FontName = Ps.FontName;
FontSize = Ps.FontSize;
FontCase = Ps.FontCase;
HShift = Ps.HShift;
VShift = Ps.VShift;

%For shifting annotations -------------------------------------------------
if Action == 2
    HiddenAnnotAx = findall(Gx, 'Tag', 'scribeOverlay');
    Txs = get(HiddenAnnotAx, 'Children');
    for j = 1:length(Txs)
        PriorUnits = get(Txs(j), 'Units');
        set(Txs(j), 'Units', 'normalized');
        TxPos = get(Txs(j), 'position');
        TxPos(1:2) = TxPos(1:2) + [HShift VShift];
        set(Txs(j), 'Position', TxPos);
        set(Txs(j), 'Units', PriorUnits);
    end
    return
end

%Removing annotations -----------------------------------------------------
HiddenAnnotAx = findall(Gx, 'Tag', 'scribeOverlay');
Txs = get(HiddenAnnotAx, 'Children');
delete(Txs);
if Action == 1; return; end %End here for deletiong only

%Adding annotations -------------------------------------------------------

%Do all operations in normalized units
Axs = findobj(Gx, 'type', 'axes');
PriorUnits = cell(length(Axs), 1);
for j = 1:length(Axs)
    PriorUnits{j} = get(Axs(j), 'Units');
    set(Axs(j), 'Units', 'normalized');
end

%Determine locations of all axes and tight insets
[UnqAxCell, AllPositions, ~] = getUniqueAxes(Gx);

%Find the columns first
GridCoord = zeros(length(UnqAxCell), 2);
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

%Determine which Axs gets what letters
LabelOrder = zeros(length(UnqAxCell), 1);
q = 0; %Label order counter adjuster
if strcmpi(Direction, 'vertical') %For each column, number
    for c = 1:max(GridCoord(:, 2))
        Idx = find(GridCoord(:, 2) == c);
        [~, SortIdx] = sort(GridCoord(Idx, 1));
        LabelOrder(Idx) = SortIdx + q;
        q = q + length(SortIdx);
    end
else
    for r = 1:max(GridCoord(:, 1))
        Idx = find(GridCoord(:, 1) == r);
        [~, SortIdx] = sort(GridCoord(Idx, 2));
        LabelOrder(Idx) = SortIdx + q;
        q = q + length(SortIdx);
    end
end

%Determine the TopLeft corners of every plot
AllOuterPositions = zeros(length(UnqAxCell), 4);
for a = 1:length(UnqAxCell)
    AllOuterPositions(a, :) = get(UnqAxCell{a}(1), 'OuterPosition');
end

%Make sure to align the labels by row and column
LabelPos = [AllOuterPositions(:, 1), sum(AllOuterPositions(:, [2 4]), 2)]; %X,Y of each label
for r = 1:max(GridCoord(:, 1))
    Idx = find(GridCoord(:, 1) == r);
    MaxY = max(LabelPos(Idx, 2));
    LabelPos(Idx, 2) = MaxY;
end
for c = 1:max(GridCoord(:, 2))
    Idx = find(GridCoord(:, 2) == c);
    MinX = min(LabelPos(Idx, 1));
    LabelPos(Idx, 1) = MinX;
end

%Get the extent via uicontrol
TempTx = uicontrol('Style', 'text');
set(TempTx, 'String', 'X)', 'FontName', FontName, 'FontSize', FontSize, 'Units', 'normalized');
TxExt = get(TempTx, 'Extent');
delete(TempTx);

%Add the annotations
for j = 1:length(UnqAxCell)
    LabelString = [char('a' + LabelOrder(j) - 1) ')'];
    if strcmpi(FontCase, 'upper')
        LabelString = upper(LabelString);
    end
    
    Tx = annotation('textbox');
    set(Tx, 'Units', 'normalized', 'LineStyle', 'none', 'FaceAlpha', 0, 'String', LabelString, 'FontSize', FontSize, 'FontName', FontName, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    TxPos = get(Tx, 'position');
    TxPos(1) = LabelPos(j, 1) - TxExt(3) + HShift;
    TxPos(2) = LabelPos(j, 2) - TxExt(4)/2 + VShift;
    TxPos(3) = TxExt(3);
    TxPos(4) = TxExt(4);
    set(Tx, 'Position', TxPos);
end

for j = 1:length(Axs)
    set(Axs(j), 'Units', PriorUnits{j});
end
