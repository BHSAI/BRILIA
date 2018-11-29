%plotHistGraph plots a default bar graph based on the inputs specified.
%
%  [Gx, Ax, Bx] =  plotHistGraph(Edges,Ncount,Normalize)
%
%  plotHistGraph(Edges,Ncount,Normalize,Param,Value,...)
%
%  plotHistGraph(Ax,Edges,Ncount,Normalize,Param,Value,...)
%
%  INPUT
%    Ax: Axes handle to plot the bar graph to. If not specified, will make
%       a new figure and axes.
%    Edges: Mx1 matrix of left edges of a histogram
%    Ncount: MxN matrix of frequency of each edge
%    Normalize ['y' or 'n']: Normalize the data
%    
%      Param           Value       Description
%      --------------  ----------  ----------------------------------------
%      Legend          Nx1 cell    Cell of legend name, in the order of the
%                                    column data of Ncount
%      Title           string      Name of the plot
%
%  OUTPUT
%    Gx: figure handle
%    Ax: axes handle
%    Bx: bar  handles

function varargout = plotHistGraph(varargin)
%Parse the initial inputs first.
if ishandle(varargin{1})
    Ax = varargin{1};
    Gx = get(Ax(1),'parent');
    varargin(1) = [];
else
    Gx = figure;
    Ax = axes;
end
Edges = varargin{1};
Ncount = varargin{2};
Normalize = varargin{3};
varargin(1:3) = [];

%Parse the varargin
P = inputParser;
addParameter(P,'Legend',{},@(x) iscell(x) || ischar(x));
addParameter(P,'Title','',@(x) ischar(x));
parse(P,varargin{:});
P = P.Results;

%Make sure edges is a Mx1 matrix
if size(Edges,1) == 1 && size(Edges,2) > size(Edges,1)
    Edges = Edges';
end

%Make sure length(Edges) == size(Ncount,1)
if length(Edges) ~= size(Ncount,1)
    if length(Edges) == size(Ncount,2)
        Ncount = Ncount';
    end
end

%Normalize the data
if lower(Normalize(1)) == 'y'
    Ncount = Ncount./repmat(sum(Ncount,1),size(Ncount,1),1);
    Yrange = [0 1];
    Yname = 'Norm. Freq.';
else
    Yrange = [0 max(Ncount(:))];
    Yname = 'Freq.';
end

%Determine the decimal units for X, and the x axis modifier
XdecRange = log10(Edges);
XdecRange(isinf(XdecRange)) = [];
Xdec = round(mean([min(XdecRange) max(XdecRange)]));
if abs(Xdec) >= 2
    Edges = Edges./10^Xdec;
    Xmodifier = sprintf('(x 10^{%d}) ',Xdec);
else
    Xmodifier = '';
end

%Determine the X range
Xincr = diff(Edges);
Xincr = Xincr(end);
Edges = [Edges; Edges(end)+Xincr];
Xrange = [min(Edges) max(Edges)];

%Determine the decimal units for Y, and the y axis modifier
YdecRange = log10(Ncount(:));
YdecRange(isinf(YdecRange)) = [];
Ydec = round(mean([min(YdecRange) max(YdecRange)]));
if abs(Ydec) >= 2 && lower(Normalize(1)) ~= 'y'
    Ncount = Ncount./10^Ydec;
    Yrange = Yrange./10^Ydec;
    Ymodifier = sprintf('(x 10^{%d}) ',Ydec);
else
    Ymodifier = '';
end

%Plot the histogram
Bx = bar(Ax,Edges(1:end-1),Ncount,'histc');
ylabel(Ax,[Ymodifier Yname]);
xlabel(Ax,Xmodifier);
set(Ax,'YLim',Yrange,'XLim',Xrange,'Xtick',Edges,'XTickLabelRotation',90); %remmeber to use 1 less tick.

%Set the legend
if ~isempty(P.Legend)
    legend(Ax,P.Legend,'Location','NorthEast')
end

%Set the title
if ~isempty(P.Title)
    title(Ax,P.Title);
end

setPlotTickDecimal(Ax,'max','max')

%Format the colors
ColorMat = [0   0   0;
            1   0   0; 
            0   1   0; 
            0.2 0.2 1;
            1   1   0; 
            0   1   1; 
            1   0   1];
for j = 1:length(Bx)
    Bx(j).FaceColor = ColorMat(j,:);
end

%Return output handles
if nargout >= 1
    varargout{1} = Gx;
    if nargout >= 2
        varargout{2} = Ax;
        if nargout >= 3
            varargout{3} = Bx;
        end
    end
end
