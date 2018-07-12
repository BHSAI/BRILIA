%labelSubplots will add the labels on the subplots at the top left corner
%of each plot. Can specify the labels to go left to right (horizontal) or
%top to bottom (vertical).
%
%  labelSubplots(Gx, 'Direction', Direction, 'FontName', FontName,
%    'FontSize', 'FontSize')
%
%  labelSubplots(Gx, 'remove')
%
%  labelSubplots(Gx, 'shift', Hshift, Vshift)
%
%  INPUT
%    Gx: figure handle 
%    Ax: axes handle that is used to find the figure handle
%    'remove': will delete the labels in Gx, if any.
%    'shift': will shift the labels in Gx, if any.
%    Hshift: shifts labels horizontally. 
%    Vshift: shifts labels vertically. 
%    
%    Param         Value            Description
%    ------------  ---------------  ---------------------------------------     
%    Direction     'horizontal'     Labels across a to z...
%                  'vertical'       Labels vertically a to z...
%    FontName      'Arial'          Font name to use for labels
%    FontSize      Integer > 0      Font size for the labels
%    FontCase      'upper'          Use upper case for labels
%                  'lower'          Use lower case for labels
%    Hshift        Double           Use for shifting labels left/right
%    Vshift        Double           Use for shifting labels up/down
%

function varargout = labelSubplots(varargin)
[Gx, varargin] = getOpenGCF(varargin{:});
P = inputParser;
addParameter(P, 'Direction', 'horizontal', @(x) ismember(lower(x), {'horizontal', 'vertical'}));
addParameter(P, 'FontName', 'Arial', @(x) ischar(x));
addParameter(P, 'FontSize', 10, @(x) isnumeric(x) && x > 2);
addParameter(P, 'FontCase', 'lower', @(x) ismember(lower(x(1:5)), {'upper', 'lower'}));
addParameter(P, 'Hshift', 0, @(x) isnumeric(x));
addParameter(P, 'Vshift', 0, @(x) isnumeric(x));

%Check if this is for removing labels
Action = 0; %0 = add labels, 1 = delete labels, 2 = shift labels
for j = 1:length(varargin)
    if ischar(varargin{j}) && ismember(lower(varargin{j}), {'remove', 'delete', 'del', 'rem'})
        Action = 1;
        varargin = {}; 
        break;
    elseif ischar(varargin{j}) && strcmpi(varargin{j}, 'shift')
        Action = 2;
        Q = inputParser;
        addRequired(Q, 'Hshift', @(x) isnumeric(x));
        addOptional(Q, 'Vshift', 0, @(x) isnumeric(x));
        parse(Q, varargin{j+1:end});
        Hshift = Q.Results.Hshift;
        Vshift = Q.Results.Vshift;
        varargin = {};
        break;
    end
end

[Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInput(P, varargin{:});
if ReturnThis
   varargout = {Ps, Pu, ExpPs, ExpPu};
   return;
end
Direction = Ps.Direction;
FontName = Ps.FontName;
FontSize = Ps.FontSize;
FontCase = Ps.FontCase;
Hshift = Ps.Hshift;
Vshift = Ps.Vshift;

if isempty(Gx)
    if isempty(findobj(0, 'type', 'figure'))
        return;
    else
        Gx = gcf;
    end
end

%For shifting annotations ---------------------------------------------------
if Action == 2
    HiddenAnnotAx = findall(Gx, 'Tag', 'scribeOverlay');
    Txs = get(HiddenAnnotAx, 'Children');
    for j = 1:length(Txs)
        PriorUnits = get(Txs(j), 'Units');
        set(Txs(j), 'Units', 'normalized');
        TxPos = get(Txs(j), 'position');
        TxPos(1:2) = TxPos(1:2) + [Hshift Vshift];
        set(Txs(j), 'Position', TxPos);
        set(Txs(j), 'Units', PriorUnits);
    end
    return;
end

%Removing annotations -----------------------------------------------------
HiddenAnnotAx = findall(Gx, 'Tag', 'scribeOverlay');
Txs = get(HiddenAnnotAx, 'Children');
delete(Txs);
if Action == 1 %End here for deletiong only
    return;
end

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
    TxPos(1) = LabelPos(j, 1) - TxExt(3) + Hshift;
    TxPos(2) = LabelPos(j, 2) - TxExt(4)/2 + Vshift;
    TxPos(3) = TxExt(3);
    TxPos(4) = TxExt(4);
    set(Tx, 'Position', TxPos);
end

%Return the axes units to prior units
for j = 1:length(Axs)
    set(Axs(j), 'Units', PriorUnits{j});
end
