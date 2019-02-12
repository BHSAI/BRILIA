%getStdColor will return a Mx3 colormap of standard colors used in figures
%for manuscripts. These are black, red, green, blue, ...

function StdColor = getStdColor(varargin)

StdColor = [0.0 0.0 0.0;
            1.0 0.0 0.0;
            0.0 0.8 0.0;
            0.0 0.0 1.0;
            0.5 0.5 0.0;
            0.5 0.0 0.5;
            0.0 0.5 0.5];

if nargin == 0; return; end
StdColor = (1 - StdColor)*(varargin{1} - 1) + StdColor;
StdColor(StdColor > 1) = 1;
StdColor(StdColor < 0) = 0;