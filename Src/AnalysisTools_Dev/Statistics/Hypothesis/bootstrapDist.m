%bootstrapDist will bootstrap a distribution based on the frequency data
%provided and provide a (1-Alpha)*100% confidence interval.
%
%  [MED, CI, MEAN, STD] = bootstrapDist(Data, varargin)
%
%  INPUT
%    Data: MxN vector of frequencies (>= 0 value) for M observations and N
%      distributions.
%
%    Param      Value      Description
%    ---------- ---------- ------------------------------------------------
%    Size       N > 0      The number of samples to draw with replacement 
%                          from each nth column distribution. Default is
%                          min(sum(Data, 1)), or the smallest count for all
%                          distributions.
%    Iter       N > 0      The number of resampling iterations to perform.
%                          Recommend > 100 to get 95% confidence intervals.
%                          Default: 100
%    Alpha      0<Alpha<1  The significance level for the conf. interval.
%    Tail      'b' 'r' 'l' The tail direction for test, b = both 2 way, 
%                          r = right tail (X>Y), and l = left tail (X<Y)
%
%  OUTPUT
%    MED:  MxN vector of median values
%    CI:   MxN cell of 1x2 double values of the (1-Alpha)*100% conf. int.
%    MEAN: MxN vector of mean values
%    STD:  MxN vector of standard deviations
%    BootData: 1xN cell of MxIter matrices of bootstrap results

function [MED, CI, MEAN, STD, BootData] = bootstrapDist(Data, varargin)
P = inputParser;
P.KeepUnmatched = 1;
addParameter(P, 'Alpha', 0.05, @(x) x >= 0 && x < 1);
addParameter(P, 'Iter', 100, @(x) x >= 0 && mod(x, 1) == 0);
addParameter(P, 'Size', [], @(x) isempty(x) || (numel(x) == 1 && x >= 0) || (numel(x) == size(Data, 2)));
addParameter(P, 'Tail', 'B', @(x) startsWith(upper(x), {'L', 'R', 'B'}));
parse(P, varargin{:});
Alpha = P.Results.Alpha;
Iter  = P.Results.Iter;
Size  = P.Results.Size;
Tail  = P.Results.Tail;

if iscell(Data)
    Data = cell2mat(Data);
end

%NOTE: for single comparison CI, L is [Alpha 1] and R is [0 1-Alpha].
%      for 2-pop comparison  CI, L is [0 1-Alpha] and R is [Alpha 1].
switch upper(Tail)
    case 'B' %2-tail both end
        Pct = [Alpha/2 1-Alpha/2]*100;
    case 'L' %Left tail
        Pct = [Alpha 1]*100;
    case 'R' %Right tail
        Pct = [0 1-Alpha]*100;
end       

if isempty(Size)
    Size = round(min(sum(Data, 1)));
end
if numel(Size) == 1
    Size = repelem(Size, 1, size(Data, 2));
end

MED = zeros(size(Data));
MEAN = zeros(size(Data));
STD = zeros(size(Data));
CI = cell(size(Data));
BootData = cell(1, size(Data, 2));
for k = 1:size(Data, 2)
    Freq = Data(:, k) / sum(Data(:, k));
    N = zeros(size(Data, 1), Iter);
    for i = 1:Iter
        A = randsample(size(Data, 1), round(Size(k)), true, Freq, varargin{:}); 
        N(:, i) = histcounts(A, [1:(size(Data, 1)+1)]);
    end

    MED(:, k) = median(N, 2);
    
    if nargout >= 2
        CIT = prctile(N, Pct, 2);
        for a = 1:size(Data, 1)
            CI{a, k} = CIT(a, :);
        end
    end

    if nargout >= 3
        MEAN(:, k) = mean(N, 2);
    end

    if nargout >= 4
        STD(:, k) = std(N, [], 2);
    end
    
    if nargout >= 5
        BootData{k} = N;
    end
end