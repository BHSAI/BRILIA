%resizeFigure is used to resize the total figure height and width, while
%scaling the internal contents. It can also be used to created an empty
%figure of specified height and dimension prior to adding plots. Always
%uses inches for the figure size.
%
%  Gx = resizeFigure(W, H)
%
%  Gx = resizeFigure(Gx, W, H)
%
%  Gx = resizeFigure(Ax, W, H)
%
%  Gx = resizeFigure(Gx, W, H, 'off', 'n')
%
%  INPUT
%    Gx: figure handle
%    Ax: axes handle that is used to determine the figure (parent of Ax)
%    H: height of the figure in inches
%    W: width of the figure in inches
%    'off' or 'on': turns Visible on or off
%    'n' or 'y': y to recenter figure after resizing, n to never recenter.
%
%  OUTPUT
%    Gx: figure handle that was modified
%
%  NOTE
%    Use resizeSubplots(Gx, ...) to resize subplots inside Gx to fill the
%    space, align axes, etc.
%
%  See also resizeSubplots
function varargout = resizeFigure(varargin)
[Gx, varargin] = getOpenGCF(varargin{:});

P = inputParser;
addParameter(P, 'FigWidth', -1, @(x) isempty(x) || isnumeric(x));
addParameter(P, 'FigHeight', -1, @(x) isempty(x) || isnumeric(x));
addParameter(P, 'Visible', '', @(x) ismember(lower(x), {'off', 'on', ''}));
addParameter(P, 'Recenter', 'y', @(x) ismember(lower(x), {'n', 'y'}));

%Decide if user used the short (x, y) or long ('Xdecimal', x, ...) input.
if ~(~isempty(varargin) && ischar(varargin{1}) && ismember(lower(varargin{1}), {'getinput', 'getinputs', 'getvarargin'}))
    IsShortInput = 1;
    for j = 1:2:length(varargin)
        if ischar(varargin{j}) && ismember(varargin{j}, {'FigWidth', 'FigHeight', 'Visible', 'getinputs', 'getinput', 'getvarargin'})
            IsShortInput = 0;
            break;
        end
    end
    if IsShortInput %Make it the longer one
        NewVarargin = {'FigWidth', -1, 'FigHeight', -1, 'Visible', 'on', 'Recenter', 'y'};
        NewVarargin(2:2:length(varargin)*2) = varargin;
        varargin = NewVarargin;
    end
end

[Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInput(P, varargin{:});
if ReturnThis
   varargout = {Ps, Pu, ExpPs, ExpPu};
   return;
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
    pause(0.05); %Concurrency bug. needs 50ms delay before changing Ax.
else %In case someone want to resize and show/hide a existing figure
    if ~isempty(Visible)
        set(Gx, 'Visible', Visible);
    end
end

%Determine valid W and H
if isempty(W) %negative means do not change
    W = -1;
end
if isempty(H) %negative means do not change
    H = -1; 
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
if W > 0
    GxPosition(3) = W;
end
if H > 0
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

if nargout >= 1
    varargout{1} = Gx;
end    
