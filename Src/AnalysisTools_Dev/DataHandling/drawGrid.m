%drawGrid will draw lines on an axes along the X and Y axis, sending the
%lines to the back as well.
%
%  Lx = drawGrid(Ax, Tick, Dir, Style, varargin)
%
%  INPUT
%    Ax: plot axes to draw on
%    Tick: tick values to add a line to
%    Dir ['x' 'y']: direction to draw the lines
%    Sytle: sytle for the line, Ex: 'k--' for black dash line
%    varargin: modifiers for the line property
%
%  OUTPUT
%    Lx: Mx1 gobject array of lines drawn
%
%  EXAMPLE
%    Ax = newplot;
%    Tick = [0:0.1:1];
%    Dir = 'x';
%    Style = 'k--';
%    Lx = drawGrid(Ax, Tick, Dir, Style, 'LineWidth', 1)

function Lx = drawGrid(Ax, Tick, Dir, Style, varargin)
Lx = gobjects(numel(Tick), 1);
hold(Ax, 'on');
ZLim = get(Ax, 'ZLim'); %Set the ticks below the ZLim.
Z = (ZLim(1) - 1) * [1 1];
switch lower(Dir)
    case 'x'
        YLim = get(Ax, 'YLim');
        for k = 1:numel(Tick)
            Lx(k) = plot(Ax, [Tick(k) Tick(k)], YLim, Style, varargin{:});
            Lx(k).ZData = Z;
        end
        set(Ax, 'YLim', YLim);
    case 'y'
        XLim = get(Ax, 'XLim');
        for k = 1:numel(Tick)
            Lx(k) = plot(Ax, XLim, [Tick(k) Tick(k)], Style, varargin{:});
        end
        set(Ax, 'XLim', XLim);
    otherwise
        error('%s: Unrecognized Dir input "%s".', mfilename, Dir);
end
hold(Ax, 'off');