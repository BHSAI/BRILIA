function mainPlotGeneUsage(varargin)
GeneUsage = getGeneUsage('', 'all');
plotGeneUsage(GeneUsage, 'maxdotsize', 300);

if nargin > 0
    savePlot(gcf, 'SaveAs', varargin{1});
else
    savePlot(gcf);
end