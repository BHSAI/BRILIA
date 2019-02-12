%plotErrorBox will plot scatter plots with error bars for multiple
%data set arranged in a MxN matrix, where N = number of group and M =
%number of categories. 
%
%  [Gx, Mx, Ex, Hx] = plotErrorBox(BarData, StdData, Param, Value, ...)
%
%  INPUT
%    BarData: MxN cell of means to plot, where col1 = data names
%
%    Param          Value
%    -------------  -------------------------------------------------------
%    'H'            MxN logical matrix where 1 = stat sig results (Can't
%                   use with 'PairH'
%    'PairH'        Mx(N*(N-1)/2) logical matrix indicating where the
%                   stat. sig. PAIR of bars exists. Note that the order of
%                   column must be according to nchoosek(1:N, K).
%    'CI'           MxN cell of 1x2 double per cell for CIs (Can't use with
%                   'STD')
%    'STD'          MxN matrix of std around the mean (Can't use with 'CI')
%    'XTickLabel'   Mx1 cell of char for names of each bar
%    'ColorOrder'   cell of colors to use in cyclic order, {'r' 'b'}
%    'IndivData'    MxZ matrix of data points for Z number of subjects.
%                   This is used to draw scatter plots instead of bars.
%    'Group'        1xZ matrix of integers for assignin each data column of
%                   IndivData into a treatment group that was ued for the
%                   BarData. Group == 1 means 1st bar of BarData, etc. 
%    'ShowIdx'      Shows only the following row idx of the input data
%                   only, to help focus on a few data points.
%
%
%  OUTPUT
%    Gx: figure handle
%    Mx: MxN cell of box objects drawn for the means
%    Ex: MxN cell of box objects drawn for the error bars
%    Hx: MxN cell of '*' objects drawn for statistical test (H) results
%  
%  NOTE:
%    Cannot use both CI and STD, so choose one.
%
%  EXAMPLE
%    BarData = [49276.9161636912 52525.491738719 20358.3837921338;2400.29417782147 5293.11687216078 4378.31911036829;905.921584647711 2588.09670290912 2365.6236482729;540.650637744852 1617.44345085827 1275.92813888296;1954.41225141057 1563.4207798055 1322.20357414472;303.544897610323 398.514689397439 303.247086859966;230.568927967006 350.666409285495 225.397766134503;192.899276854432 160.872131123449 161.470072140876;315.798134791448 149.719231276655 211.190167611722;116.278247095945 35.0055469550298 86.7144462447098;115.39717142568 37.9562558888552 64.7596119309166];
%    Stds = [2508.44143062972 1333.03977074255 6841.65563332658;225.186939081426 851.178176130221 1848.43383686642;157.789679430599 528.173776068646 1022.64592116386;86.2941084743158 265.320373082552 469.521776950312;171.803873310362 368.344978299071 311.974058069722;48.8708802492588 106.040467505818 60.9289040997831;34.5558985641037 85.8861423746409 57.2349187132005;26.5797031375654 40.9543226306482 78.8101304930925;36.3204221471561 11.80916160975 104.151516332017;23.3593669210256 23.9285917291214 65.3515853647744;19.2135607533072 14.5652242219349 56.334262924102];
%    H = [0 0 1;0 1 0;0 1 0;0 1 0;0 0 0;0 0 0;0 0 0;0 0 0;0 1 0;0 1 0;0 1 0];
%
%    [Gx, Mx, Ex, Hx] = plotErrorBox(BarData, 'STD', Stds, 'H', H);


function [Gx, Mx, Ex, Hx] = plotErrorBox(varargin)
%Initial parsing of axes
if isempty(varargin)
    error('%s: Not enough inputs.', mfilename);
end
if nargin >= 1
    if isa(varargin{1}, 'matlab.graphics.axis.Axes')
        if isvalid(varargin{1})
            Ax = varargin{1};
            Gx = get(Ax, 'parent');
        end
        varargin = varargin(2:end);
    elseif isa(varargin{1}, 'matlab.ui.Figure')
        if isvalid(varargin{1})
            Gx = varargin{1};
            Ax = findobj(Gx, 'type', 'axes');
            if isempty(Ax)
                Ax = axes(Gx);
            else
                Ax = Ax(1);
            end
        end
        varargin = varargin(2:end);
    elseif isempty(varargin{1})
        Gx = figure;
        Ax = axes(Gx);
        varargin = varargin(2:end);
    else
        Gx = figure; 
        Ax = axes(Gx);
    end
end    
%     if isempty(varargin{1}) || (~ischar(varargin{1}) && ~isnumeric(varargin{1}) && ~ishandle(varargin{1}))
%         Gx = figure;
%         Ax = axes(Gx);
%         varargin = varargin(2:end);
%     elseif isa(varargin{1}, 'matlab.graphics.axis.Axes')
%         Ax = varargin{1};
%         Gx = get(Ax, 'parent');
%         varargin = varargin(2:end);
%     elseif isa(varargin{1}, 'matlab.ui.Figure')
%         if ~isvalid(varargin{1}) %deleted figure handle. recreate.
%             Gx = figure;
%             Ax = axes(Gx);
%         else
%             Gx = varargin{1};
%             Ax = get(Gx, 'children');
%         end
%         varargin = varargin(2:end);
%     else
%         Gx = figure;
%         Ax = axes(Gx);
%     end
% else
%     Gx = figure;
%     Ax = axes(Gx);
% end
BarData = varargin{1};
varargin = varargin(2:end);
[M, N] = size(BarData);

%Options parsing
P = inputParser;
addParameter(P, 'H',  [], @(x) ismatrix(x) && isequal([M, N], size(x)));
addParameter(P, 'PairH',  [], @(x) isempty(x) || (ismatrix(x) && isequal(M, size(x, 1)))); % && size(x, 2) == N*(N-1)/2));
addParameter(P, 'STD',[], @(x) ismatrix(x) && isequal(size(BarData), size(x)));
addParameter(P, 'CI', [], @(x) iscell(x) && all(cellfun(@(b) numel(b) == 2, x(:))));
addParameter(P, 'SigMarkerOrder', {'x' 'o' '+' '*' 's' 'd' 'p' 'h' '<' '^' '>'}, @(x) iscell(x) && all(cellfun(@ischar, x)));
addParameter(P, 'MarkerOrder', {'o' '^' 's' 'd' 'p' 'h' '<' '>'}, @(x) iscell(x) && all(cellfun(@ischar, x)));
addParameter(P, 'ColorOrder',  {[0.2 0.2 0.2] [0.8 0.2 0.2] [0.2 0.8 0.2] [0.2 0.2 0.8] 'm' 'c'}, @(x) iscell(x) && all(cellfun(@ischar, x)));
addParameter(P, 'XTickLabel', 1:size(BarData, 1), @(x) (ismatrix(x) || iscell(x)) && numel(x) == size(BarData, 1));
addParameter(P, 'IndivData', [], @(x) ismatrix(x) && size(BarData, 1) == size(x, 1));
addParameter(P, 'Group', [], @(x) ismatrix(x));
addParameter(P, 'MaxW', 0.8, @(x) isnumeric(x) && x > 0 && x <= 1);
addParameter(P, 'ShowIdx', [], @isnumeric);
parse(P, varargin{:});
H       = P.Results.H;
PairH   = P.Results.PairH;
STD     = P.Results.STD;
CI      = P.Results.CI;
SigMkr  = P.Results.SigMarkerOrder; %SigMkrNum = 1;
Mkr     = P.Results.MarkerOrder; %MkrNum = 1;
Clr     = P.Results.ColorOrder; %ClrNum = 1;
XTickLabel = P.Results.XTickLabel;
IndivData  = P.Results.IndivData;
Group   = P.Results.Group;
MaxW    = P.Results.MaxW;
ShowIdx = P.Results.ShowIdx;

%Filter data first
if ~isempty(ShowIdx)
    if ~isempty(BarData)
        BarData = BarData(ShowIdx, :);
    end
    if ~isempty(CI)
        CI = CI(ShowIdx, :);
    end
    if ~isempty(STD)
        STD = STD(ShowIdx, :);
    end
    if ~isempty(IndivData)
        IndivData = IndivData(ShowIdx, :);
    end
    if ~isempty(H)
        H = H(ShowIdx, :);
    end
    if ~isempty(PairH)
        PairH = PairH(ShowIdx, :);
    end
    if ~isempty(XTickLabel)
        XTickLabel = XTickLabel(ShowIdx, :);
    end
end

%Compute the left edges of the rectangles
RectW = MaxW/size(BarData, 2);
Xshift = (0:RectW:(size(BarData, 2)-1)*RectW);
Xshift = Xshift - mean(Xshift);
X = [1:size(BarData, 1)]'; %#ok<NBRAK>
RectL = X + Xshift - RectW/2;

%Compute CI related stuff
if ~isempty(STD) && ~isempty(CI)
    error('%s: Cannot specify both a STD and a CI. Choose one.', mfilename);
elseif isempty(STD) && isempty(CI)
    DrawErrorBox = 0;
else
    DrawErrorBox = 1;
    if ~isempty(STD)
        CI = cell(size(BarData));
        for a = 1:numel(BarData)
            CI{a} = [BarData(a)-STD(a) BarData(a)+STD(a)];
        end
    end
    RectH = cellfun(@(x) diff(x), CI);
    RectB = cellfun(@(x) x(1), CI);
end

%Draw CI error
Ex = gobjects(size(BarData));
if DrawErrorBox
    ClrNum = 1;
    for a = 1:numel(BarData)
        Ex(a) = rectangle('position', [RectL(a), RectB(a), RectW, RectH(a)], 'FaceColor', 'none', 'EdgeColor', 'none');
    end
    for c = 1:size(BarData, 2)
        if isempty(IndivData)
            set(Ex(:, c), 'FaceColor', Clr{ClrNum})
        else
            set(Ex(:, c), 'EdgeColor', Clr{ClrNum})
        end
        ClrNum = incr(ClrNum, length(Clr), 1);
    end
end

%Draw the BarData IF indiv data is not specified
Mx = gobjects(size(BarData));
if isempty(IndivData)
    ClrNum = 1;
    for a = 1:numel(BarData)
        Mx(a) = rectangle(Ax, 'position', [RectL(a), 0, RectW, BarData(a)], 'FaceColor', 'none', 'EdgeColor', 'k');
    end
    if ~DrawErrorBox %If there's no STD/CI, color code bar instead
        for c = 1:size(BarData, 2)
            set(Mx(:, c), 'FaceColor', Clr{ClrNum});
            ClrNum = incr(ClrNum, length(Clr), 1);
        end
    end
end

if ~isempty(IndivData) && ~isempty(Group)
    UnqGroup = unique(Group);
    hold(gca, 'on')
    X = 1:size(IndivData, 1);
    Sx = gobjects(length(UnqGroup), 1);
    MkrNum = 1;
    ClrNum = 1;
    for y = 1:length(UnqGroup)
        GrpLoc = UnqGroup(y) == Group;
        Yall = IndivData(:, GrpLoc);
        Yall = Yall(:);
        Xall = repmat(X(:), 1, sum(GrpLoc));
        Xall = Xall(:);
        
        Sx(y) = scatter(Ax, Xall+Xshift(y) , Yall, 20, Clr{ClrNum}, Mkr{MkrNum}, 'fill');
        MkrNum = incr(MkrNum, length(Mkr), 1);
        ClrNum = incr(ClrNum, length(Clr), 1);
    end
    hold(gca, 'off');
end

%Draw a * above the significant ones
Hx = gobjects(size(BarData));
if ~isempty(H)
    ClrNum = 1;
    for a = 1:numel(BarData)
        if ~isnan(H(a)) && H(a)
            Ypos = 0.02*max(max(RectH + RectB)) + (RectB(a) + RectH(a));
            Hx(a) = text(Ax, RectL(a)+RectW/2, Ypos, '*', 'FontName', 'Courier New', 'FontSize', 14, 'HorizontalAlignment', 'center');
        end
    end
    for c = 1:size(BarData, 2)
        set(Hx(H(:, c)>0, c), 'Color', Clr{ClrNum});
        ClrNum = incr(ClrNum, length(Clr), 1);
    end
end

%Determine how much to shift the above-bar labels
if ~isempty(IndivData)
    MaxY = max(max(IndivData, [], 2));
elseif ~isempty(BarData)
    MaxY = max(max(BarData, [], 2));
else
    YLim = get(gca, 'YLim');
    MaxY = YLim(2);
end
Yshift = 0.01*MaxY;

%Draw the pairwise difference
Px = gobjects(size(BarData, 1), 1); %Stores the text
Lx = gobjects(size(BarData, 1), 1); %Stores the divider lines
MaxYLim = 0;
if ~isempty(PairH)
    Combo = nchoosek(1:size(BarData, 2), 2);
    if DrawErrorBox
        Y = max(RectB + RectH, [], 2);
    else
        Y = max(BarData, [], 2);
    end
    TxtMkr = cell(1, size(PairH, 2));
    SigMkrNum = 1;
    for k = 1:size(PairH, 2)
        TxtMkr{k} = SigMkr{SigMkrNum};
        SigMkrNum = incr(SigMkrNum, length(SigMkr), 1);
    end
    hold(gca, 'on')
    for r = 1:size(PairH, 1)
        NumSig = sum(PairH(r, :));
        if NumSig > 0
            SigIdx = find(PairH(r, :) > 0);
            DrawMkr = cell(1, length(SigIdx));
            for z = 1:length(SigIdx)
                DrawMkr{z} = sprintf('%dv%d', Combo(SigIdx(z), 1), Combo(SigIdx(z), 2));
            end
            Ycoor = Y(r)+Yshift*2;
            Px(r) = text(Ax, r, Ycoor, DrawMkr, 'FontName', 'Courier New', 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            Lx(r) = plot(Ax, [r+Xshift(1)-RectW/2 r+Xshift(end)+RectW/2], Yshift+Y(r)*[1 1], 'k-', 'LineWidth', 1);
            if Ycoor > MaxYLim
                MaxYLim = Ycoor;
            end
        end
    end
    hold(gca, 'off');
end

%Final formatting of figures
XLim = [1-MaxW size(BarData, 1)+MaxW];
YLim = get(gca, 'YLim');
YTick = get(gca, 'YTick');
DiffYTick = diff(YTick(1:2));
if YLim(2)-DiffYTick < MaxYLim
    YLim(2) = YLim(2) + 2*DiffYTick;
end
ylabel('Frequency');
set(Ax, 'XLim', XLim, 'YLim', YLim, 'FontName', 'Courier New', 'FontSize', 10, 'XTick', 1:size(BarData, 1));
if ~isempty(XTickLabel)
    set(Ax, 'XTickLabel', XTickLabel, 'XTickLabelRotation', 90);
end
resizeSubplots(gca);