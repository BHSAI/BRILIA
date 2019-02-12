%calcPairChiStat will compute the chi sq and p-values of distribution
%matches between frequency data sets, arranged in a M x N matrix, where M
%rows are for a bin, and N columns are for the number of distributions.
%
%  [ChiMat, PMat] = calcPairChiStat(FreqMat)
%
%  INPUT
%    FreqMat: M (bins) x N (frequency data) matrix
%
%  OUTPUT
%    ChiMat: NxN diagonally symmetric matrix of Chi-Square statistics value
%    PMat: NxN diagonally symmetric matrix of p-values

function [ChiMat, PMat] = calcPairChiStat(FreqMat)
ChiMat = zeros(size(FreqMat, 2));
PMat   = zeros(size(FreqMat, 2));
for r = 1:size(FreqMat, 2)
    for c = 1:size(FreqMat, 2)
        Freq2 = rescaleDist(FreqMat(:, r), FreqMat(:, c));
        [ChiMat(r, c), PMat(r, c)] = calcChiStat(FreqMat(:, r), Freq2, 'norm');
    end
end
