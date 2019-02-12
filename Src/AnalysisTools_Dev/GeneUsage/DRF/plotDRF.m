%plotDRF will plot the D reading frame frequency data stored in the data
%structure S.
%
%  plotDRF(S, SizeFilter)
%
%  plotDRF(S, SizeFilter, 'norm')
%
%  INPUT
%    S: structure of VDJdata and VDJheader
%    SizeFilter ['AC' 'BC' 'BCN' 'TOPN' 'BOTN']: filter for the group based
%         on clonotype sizes
%      AC - all clonotypes
%      BC - branched clonotypes with >= 2 unique sequences per clonotype
%      BCN - branched clonotypes with >= N unique sequences per clonotype
%      TOPN - top N clonotypes based on total Template count per clonotype
%      BOTN - top N clonotpyes based on total Template count per clonotype
%   'norm': will normalize the plots

function plotDRF(S, SizeFilter, varargin)
if ~isempty(varargin) && any(startsWith(varargin, 'norm', 'ignorecase', true))
    Normalize = true;
else
    Normalize = false;
end

ValidFiltType = {'BC', 'AC', 'SC', 'TOP', 'BOT'};
if ~startsWith(SizeFilter, ValidFiltType)
    error('%s: Invalid option for SizeFilter, "%s".', mfilename, SizeFilter);
end

D = cell(1, length(S));
for f = 1:length(S)
    D{f} = formatDRF(S(f).([SizeFilter '_DRF']), 4);
end
C = conformDist(D{:});

XLabel = C(:, 1);
for f = 1:length(S)
    Y = cell2mat(C(:, 1+[(f-1)*3+1:f*3]));
    if Normalize
        Y = Y./repmat(sum(Y, 2), 1, size(Y, 2));
    end
    figure
    Bx = bar(Y, 'stacked');
    Bx(1).FaceColor = 'r';
    Bx(2).FaceColor = 'g';
    Bx(3).FaceColor = 'b';
    set(gca, 'XTick', 1:length(XLabel), 'XTickLabel', XLabel, 'XTickLabelRotation', 90, 'XLim', [0 length(XLabel)+1])
    ylabel('Fraction');
    title(getTitleFileName(S(f).FileName));
    resizeFigure(gcf, 5, 4, 'on', 'n');
    resizeSubplots(gcf)
    
    [~, FileTemp] = parseFileName(S(f).FileName);
    savePlot(gcf, 'SaveAs', [FileTemp(1:end-4) '.DRF.png']);
end