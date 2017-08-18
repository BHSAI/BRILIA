%centerFigureOnMonitor will center figures within the monitor. This is
%used mainly to prevent code-generated figures from being shown off the
%screen, making it hard for users to get access to it.
%
%  centerFigureOnMonitor(Gx)
%
%  centerFigureOnMonitor(Ax)
%
%  centerFigureOnMonitor
%
%  INPUT
%    Gx: figure handle to center. If empty, operates on all open figures.
%        Gx = get(0, 'children')
%    Ax: axes handle used to get the figure handle to center
function centerFigureOnMonitor(varargin)
Gxs = getOpenGCF(varargin{:});
if isempty(Gxs)
    return;
end

for g = 1:length(Gxs)
    Gx = Gxs(g);

    %Reposition the figure at the center first.
    OrigGxUnits = get(Gx, 'Units');
    set(Gx, 'Units', 'pixels');
    FigureOuterPos = get(Gx, 'outerposition'); %Includes the menu, etc
    ScreenSize = get(0, 'screensize');

    NewFigureOuterPos = FigureOuterPos;
    NewFigureOuterPos(1) = ScreenSize(3)/2 - FigureOuterPos(3)/2;
    NewFigureOuterPos(2) = ScreenSize(4)/2 - FigureOuterPos(4)/2;

    %If the figure top is going off the screen, move down.
    CurTop = sum(NewFigureOuterPos([2 4]));
    if CurTop > ScreenSize(4)
        Adj = ScreenSize(4) - CurTop; 
        NewFigureOuterPos(2) = NewFigureOuterPos(2) + Adj;
    end

    set(Gx, 'outerposition', NewFigureOuterPos);
    set(Gx, 'Units', OrigGxUnits);
end
