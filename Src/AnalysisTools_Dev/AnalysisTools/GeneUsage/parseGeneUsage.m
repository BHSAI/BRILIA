%parseGeneUsage will add the frequency numbers from the N-D matrix from
%getGeneUsage to make it easier to plot frequency distribution based on a
%V, D, or J gene.
%
%  [FreqData, GeneName] = parseGeneUsage(GetUsage, Gene)
%
%  INPUT
%    GeneUsage: structure results from the getGeneUsage.m 
%    Gene ['V', 'D', 'J']: the gene frequency to return
%
%  OUTPUT
%    FreqData: Nx1 frequency 
%    GeneName: the gene names
%
%  EXAMPLE
%    GeneUsage = getGeneUsage('File1.csv', 'all')
%    [FreqData, GeneName] = parseGeneUsage(GeneUsage, 'V')
%    barh(FreqData)
%    GeneName = strrep(GeneName, 'IGH', '')
%    set(gca, 'YTick', 1:length(GeneName), 'YTickLabel', GeneName, 'YLim', [0 length(GeneName)+1])
%    resizeSubplots(gca)
function [FreqData, GeneName] = parseGeneUsage(GeneUsage, Gene)

switch upper(Gene(1))
    case 'V'
        FreqData = sum(sum(GeneUsage.Count, 2), 3);
        GeneName = GeneUsage.Names{1};
    case 'D'
        FreqData = sum(sum(GeneUsage.Count, 1), 3);
        GeneName = GeneUsage.Names{2};
    case 'J'
        FreqData = sum(sum(GeneUsage.Count, 1), 2);
        GeneName = GeneUsage.Names{3};
end
FreqData = FreqData(:, 1);