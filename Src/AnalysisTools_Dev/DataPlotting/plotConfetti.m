%plotConfetti will draw the confetti image given either a list of group
%sizes or the 2D confetti map matrix.
%
%  INPUT
%    GrpSizes: Nx1 array of sizes of groups
%    Confetti: sqrt(N) x sqrt(N) area matrix of unique group numbers
%              (from makeConfettiMap.m)
%
%  OUTPUT
%    Gx: figure handle of the plot
%
%  EXAMPLE
%    Sizes = [100 30 20 10 20 1 1 3 2 4];
%    plotConfetti(Sizes)

function [Gx, Confetti] = plotConfetti(varargin)
%Using persistent rand to prevent generating a different figure every time
persistent s
if isempty(s)
    s = RandStream('mt19937ar', 'Seed', 1);
else
    s.reset;
end

[Axs, varargin] = getOpenGCA(varargin{:});
if isempty(varargin)
    error('%s: Incorrect input. Need a matrix of sizes.', mfilename);
else
    Confetti = varargin{1};
    if numel(varargin) == 2
        Color = varargin{2};
    else
        Color = [];
    end
end
if min(size(Confetti)) == 1 && ~(size(Confetti, 1) == size(Confetti, 2)) %You're NOT given confetti map (but rather the Sizes). Compute confetti map.
    Confetti = makeConfettiMap(Confetti);
end
if isempty(Axs)
    Gx = figure;
    Axs = gca;
end
image(Axs(1), Confetti, 'CDataMapping', 'direct');
if isempty(Color)
    colormap(Axs(1), rand(s, max(Confetti(:)), 3));
else
    colormap(Axs(1), Color);
end
setAxes(Axs(1), 'XTickLabel', '', 'YTickLabel', '');

%Make sure Axs(1) is NOT a subplot of a bigger picture
OtherAxs = get(get(Axs(1), 'parent'), 'children');
if numel(OtherAxs) == 1
    resizeFigure(Axs(1), 5, 5);
    resizeSubplots(Axs(1));
    centerFigureOnMonitor(Axs(1))
end






