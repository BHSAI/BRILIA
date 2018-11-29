%getMaxPlotPosition will get the maximum plot position possible without
%running the x and y labels off the defined OuterPos area.
%
%  MaxPos = getMaxPlotPosition(Ax, OuterPos)
%
%  MaxPos = getMaxPlotPosition(Ax, OuterPos, Param, Value)
%
%  INPUT
%    Ax: axis handle
%    OuterPos: 1x4 matrix of the axes outer position you want Ax to fit
%      within, given as [Xloc, Yloc, Width, Height]
%    Param        Value             Description
%    ------------ ----------------- ---------------------------------------
%    LeftPad      Double {0, 1}     Add left border, or preserve (-1)
%    RightPad     Double {0, 1}     Add right border, or presreve (-1)
%    TopPad       Double {0, 1}     Add top border, or preserve (-1)
%    BottomPad    Double {0, 1}     Add bottom border, or preserve (-1)
%
%  OUTPUT
%    MaxPos: 1x4 matrix of the axes position that will fill the entire
%      region specified by OuterPos
%
%  NOTE
%    All resizing uses the normalized axes units to fill the figure space.
%    All axes units are restored to the original units after resizing.
%
function varargout = getMaxPlotPosition(Ax, OuterPos, varargin)
if ishandle(Ax) && ~strcmpi(class(Ax), 'matlab.graphics.axis.Axes')
    error('%s: Need to input a valid axes handle.', mfilename);
end
P = inputParser;
addParameter(P, 'LeftPad', 0, @(x) isnumeric(x));
addParameter(P, 'RightPad', 0, @(x) isnumeric(x));
addParameter(P, 'TopPad', 0, @(x) isnumeric(x));
addParameter(P, 'BottomPad', 0, @(x) isnumeric(x));

[Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInput(P, varargin{:});
if ReturnThis
   varargout = {Ps, Pu, ExpPs, ExpPu};
   return;
end
PadL = Ps.LeftPad;
PadR = Ps.RightPad;
PadT = Ps.TopPad; 
PadB = Ps.BottomPad;

%Calculate the tight inset
TightPos = get(Ax, 'TightInset');
MaxPos = zeros(1, 4);
MaxPos(1) = OuterPos(1) + TightPos(1) + PadL;
MaxPos(2) = OuterPos(2) + TightPos(2) + PadB;
MaxPos(3) = OuterPos(3) - TightPos(1) - TightPos(3) - PadL - PadR;
MaxPos(4) = OuterPos(4) - TightPos(2) - TightPos(4) - PadT - PadB;

if nargout >= 1
    varargout{1} = MaxPos;
end
