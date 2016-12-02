%plotGeneMods will take a cell array of numbers Ndata, and plot the
%histograms.
%
%plotGeneMods(Ndata,Edges,Headers) will force and header
%returns the plot bar handle
function [H1, Ax] = plotGeneMods(Ndata,varargin)
%Determine the edges
if isempty(varargin) || (length(varargin)>=1 && isempty(varargin{1})) %Autofind edges
    MaxNum = 0;
    MinNum = 0;
    for j = 1:length(Ndata)
        MaxNumNow = max(Ndata{j});
        if MaxNumNow > MaxNum
            MaxNum = MaxNumNow;
        end
        MinNumNow = max(Ndata{j});
        if MinNumNow < MinNum
            MinNum = MinNumNow;
        end
    end
    Edges = MinNum:MaxNum;
else
    Edges = varargin{1};
end

%Determine the legend
Headers = [];
if length(varargin) >= 2
    Headers = varargin{2};
end
if isempty(Headers)
    Headers = cell(1,length(Ndata));
    for k = 1:length(Ndata)
        Headers{k} = sprintf('Data %d',k);
    end
end

Normalize = 0;
if length(varargin) == 3
    if strcmpi(varargin{3},'norm');
        Normalize = 1;
    end
end

%Calculate frequency histogram
Ncount = zeros(length(Edges),length(Ndata));
for j = 1:length(Ndata)
    Ncount(1:end-1,j) = histcounts(Ndata{j},Edges)';
end

if Normalize == 1
    Ncount = Ncount./repmat(sum(Ncount,1),size(Ncount,1),1);
    Yrange = [0 1];
else
    Yrange = [0 max(Ncount(:))];
end

%Plot the histogram
figure
H1 = bar(Edges,Ncount);
Ax = gca;
legend(Headers);
ylim(Ax,Yrange)
xlim(Ax,[Edges(1)-1  Edges(end)+1]);
set(Ax,'XTickLabel',Edges);
set(Ax,'XTick',Edges);
set(Ax,'XTickLabelRotation',70);
set(Ax,'OuterPosition',[0 0 1 1]);

%Format the colors
ColorMat = [1 0 0; 0 1 0; 0.2 0.2 1;...
            1 1 0; 0 1 1; 1 0 1];
for j = 1:length(H1)
    H1(j).FaceColor = ColorMat(j,:);
end

xlabel(Ax,'Counts')
if Normalize == 1
    ylabel(Ax,'Norm Freq');
else
    ylabel(Ax,'Freq');
end


