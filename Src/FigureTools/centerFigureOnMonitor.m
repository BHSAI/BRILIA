%centerFigureOnMonitor will center figures on the screen. This is
%used mainly to prevent code-generated figures from being shown off the
%screen, making it hard for users to get access to it.
%
%  centerFigureOnMonitor
%
%  centerFigureOnMonitor(Gx)
%
%  INPUT
%    Gx: figure or axes handle. If empty, operates on all open figures.
%
%  EXAMPLE
%    Gx = loadDemoPlot;
%    centerFigureOnMonitor(Gx);
%
function centerFigureOnMonitor(varargin)
Gxs = getOpenGCF(varargin{:});
if isempty(Gxs); return; end

ScreenSize = get(0, 'screensize');
ScreenWCen = ScreenSize(3) / 2;
ScreenHCen = ScreenSize(4) / 2;
for g = 1:length(Gxs)
    Gx = Gxs(g);
    OrigGxUnits = Gx.Units;
    Gx.Units = 'pixels';

    FigWCen = round(Gx.OuterPosition(3) / 2);
    FigHCen = round(Gx.OuterPosition(4) / 2);

    OverTop = max(0, FigHCen + ScreenHCen - ScreenSize(4));
    NewOuterPosition = Gx.OuterPosition;
    NewOuterPosition(1) = ScreenWCen - FigWCen;
    NewOuterPosition(2) = ScreenHCen - FigHCen - OverTop;   

    Gx.OuterPosition = NewOuterPosition;
    Gx.Units = OrigGxUnits;
end