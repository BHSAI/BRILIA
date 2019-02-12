%correctP will use multiple-hypothesis testing correction to identify
%statistically signficant p-values.
%
%  Idx = correctP(P, Alpha, Method)
%
%  INPUT
%    P: Nx1 vector of p-values for N tests
%    Alpha: signficance level you want
%    Method: method for multiple-hypothesis testing. Default is 'b'.
%      'b'  - Bonferonni 
%             1) Find all P that satisfies P < (Alpha/N)
%      'bh' - Benjamini-Hochberg 
%             1) Sort P and find last P, Pmax, that satisfies 
%                P < (i/N)*Alpha) , where i = rank #
%             2) Find all P less than Pmax
%      'mbh' - Modified Benjamini-Hochberg by Rom 1990, explained by Rand
%              R. Wilcox 2012 paper
%
%  OUTPUT
%    Idx: Index of test numbers that have p < p_corrected
%
%  EXAMPLE
%    P = [0.001 0.039 0.008 0.041 0.042 0.060 0.07];
%    Alpha = 0.05;
%    Idx = correctP(P, Alpha, 'bh')
%    Idx =
%         1
%         3
function Idx = correctP(P, Alpha, Method)
if nargin < 2
    Alpha = 0.05;
    if nargin < 3
        Method = 'b';
    end
end

switch lower(Method)
    case 'b' %Bonferroni correction
        Idx = find(P <= (Alpha/length(P)));
    case 'bh' %Benjamini-Hochberg
        [Prank, Pidx] = sort(P(:));
        Ptest = Alpha * [1:numel(P)]'/numel(P);
        PPidx = find(Prank <= Ptest, 1, 'last');
        Idx = Pidx(1:PPidx);
    case 'mbh' %Modified Benjamini-Hochberg by Rom 1990, explained by Rand R. Wilcox 2012 paper. The sequential rejection method.
        [Prank, Pidx] = sort(P(:), 'descend');
        k = 1;
        for j = 1:length(Prank)
            if Prank(k) <= Alpha/k
                Idx = Pidx(k:end);
                break
            end
            k = incr(k);
        end
    otherwise
        error('%s: Unrecognized Method "%s".', mfilename, Method);
end