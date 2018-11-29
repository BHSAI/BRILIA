%resizeFigure is used to resize the total figure height and width, while
%scaling the internal contents. It can also be used to created an empty
%figure of specified height and dimension prior to adding plots. Always
%uses inches for the figure size.
%
%  Gx = resizeFigure(W, H)
%
%  Gx = resizeFigure(Gx, W, H)
%
%  Gx = resizeFigure(Gx, W, H, Visible, Recenter)
%
%  Gx = resizeFigure(Gx, Param, Value, ...)
%
%  INPUT
%    Gx: figure or axes handle handle. If empty, will use current figure
%    H: height of the figure in inches
%    W: width of the figure in inches
%    Visible ['off' 'on']: turns Visible on or off
%    Recenter ['n' 'y']: yes or no to recenter figure
%
%    Parameter    Values(*)   Description
%    -----------  ----------  ---------------------------------------------
%    FigWidth     * -1        Do not modify figure width
%                   double    Set figure width in inch
%    FigHeight    * -1        Figure height in inch. Negative = do not edit
%                   double    Set figure height in inch
%    Visible      * 'on'      Show the figure
%                   'off'     Hide the figure
%    Recenter     * 'n'       Do not recenter figure on screen after resize
%                   'y'       Recenter figure on screen after resize
%
%  OUTPUT
%    Gx: figure handle that was modified
%
%  EXAMPLE
%    Gx = loadDemoPlot;
%    resizeFigure(Gx, 4, 4);
%
function varargout = resizeFigure(varargin)
[Gx, varargin] = getOpenGCF(varargin{:});

P = inputParser;
addParameter(P, 'FigWidth', -1, @(x) isempty(x) || isnumeric(x));
addParameter(P, 'FigHeight', -1, @(x) isempty(x) || isnumeric(x));
addParameter(P, 'Visible', '', @(x) ismember(lower(x), {'off', 'on', ''})); %empty as in leave unchanged
addParameter(P, 'Recenter', 'n', @(x) ismember(lower(x), {'n', 'y'}));
if ~any(cellfun('isclass', varargin, 'char')) %Using short input style
    NewVarargin = {'FigWidth', -1, 'FigHeight', -1, 'Visible', '', 'Recenter', 'n'};
    NewVarargin(2:2:length(varargin)*2) = varargin;
    varargin = NewVarargin;
end
[Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInput(P, varargin{:});
if ReturnThis
   varargout = {Ps, Pu, ExpPs, ExpPu};
   return
end
W = Ps.FigWidth;
H = Ps.FigHeight;
Visible = Ps.Visible;
Recenter = Ps.Recenter;

if isempty(Gx)
    if isempty(Visible)
        Visible = 'on';
    end
    Gx = figure('Visible', Visible);
    axes(Gx);
    drawnow('nocallbacks');
else %In case someone want to resize and show/hide a existing figure
    if ~isempty(Visible)
        set(Gx, 'Visible', Visible);
    end
end

%Convert all axes to normal, saving original units
Axs = findobj(Gx, 'type', 'axes');
OrigAxUnits = cell(length(Axs), 1);
for j = 1:length(Axs)
    OrigAxUnits{j} = get(Axs(j), 'Units');
    set(Axs(j), 'Units', 'normalized');
end

%Resize the figure now.
OrigGxUnits = get(Gx, 'Units');
set(Gx, 'Units', 'inches');
GxPosition = get(Gx, 'Position');
if ~isempty(W) && W > 0
    GxPosition(3) = W; 
end
if ~isempty(H) && H > 0
    GxPosition(4) = H; 
end
set(Gx, 'Position', GxPosition);
set(Gx, 'PaperPosition', [0 0 GxPosition(3:4)]);
set(Gx, 'Units', OrigGxUnits);

%Reset axes units to orignal units
for j = 1:length(Axs)
    set(Axs(j), 'Units', OrigAxUnits{j});
end

if strcmpi(Recenter, 'y')
    centerFigureOnMonitor(Gx);
end

varargout{1} = Gx;