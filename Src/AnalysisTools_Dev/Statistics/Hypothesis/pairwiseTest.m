%pairwiseTest will test for difference of means based on 2 vectors of
%values, returning a CI and P value for the observed difference in means.
%
%  [H, P, CI] = pairwiseTest(X, Y, varargin)
%
%  INPUT
%    X: vector of values from sample A
%    Y: vector of values from sample B
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
%    Alpha      0 < a < 1  The significance level for the conf. interval.
%    Tail      'b' 'r' 'l' The tail direction for test, b = both 2 way, 
%                          r = right tail (X>Y), and l = left tail (X<Y)
%
%  OUTPUT
%    H: 1 for reject H0, 0 for fail to reject H0
%    P: p-value for observing difference in means of X and Y
%    CI: 1x2 matrix for the (1-Alpha)*100% confidence interval of the
%      difference in medians
%    DXY: all combo difference between X and Y (X-Y)
%  EXAMPLE
%    X = randn(1, 100)+5;
%    Y = randn(1, 50)+4;
%    [H, P, CI] = ttest2(X, Y, 'alpha', 0.05);
%    [H, P, CI] = pairwiseTest(X, Y, 'alpha', 0.05);

function [H, Pv, CI, DXY] = pairwiseTest(X, Y, varargin)
P = inputParser;
addParameter(P, 'Alpha', 0.05, @(x) x >= 0 && x < 1);
addParameter(P, 'Tail', 'B', @(x) startsWith(upper(x), {'L', 'R', 'B'}));
addParameter(P, 'CombVec', [], @isnumeric);
parse(P, varargin{:});
Alpha   = P.Results.Alpha;
Tail    = P.Results.Tail;
CombVec = P.Results.CombVec;

%NOTE: for single comparison CI, L is [Alpha 1] and R is [0 1-Alpha].
%      for 2-pop comparison  CI, L is [0 1-Alpha] and R is [Alpha 1].
switch upper(Tail)
    case 'B' %2-tail both end
        Pct = [Alpha/2 1-Alpha/2]*100;
    case 'L' %Left tail
        Pct = [0 1-Alpha]*100;
    case 'R' %Right tail
        Pct = [Alpha 1]*100;
end      

if isempty(CombVec)
    CombVec = combvec2(1:length(X), 1:length(Y));
end
if size(CombVec, 1) > size(CombVec, 2)
    DXY = X(CombVec(:, 1)) - Y(CombVec(:, 2));
else
    DXY = X(CombVec(1, :)) - Y(CombVec(2, :));
end
DXY = sort(DXY);

MedianDXY = median(DXY);
switch upper(Tail)
    case 'B'%Any difference
        if median(DXY) == 0 %If  you're right on the middle, P = 1
            Pv = 1;
        else
            if MedianDXY > 0
                TailNum = sum(DXY < 0) + 0.5*sum(DXY == 0);
            else
                TailNum = sum(DXY > 0) + 0.5*sum(DXY == 0);
            end
            if TailNum == 0 %just set to smallest value
                Pv = realmin;
            else
                Pv = 2*TailNum / numel(DXY);
            end
        end
    case 'L' %(X < Y)
        TailNum = sum(DXY > 0) + 0.5*sum(DXY == 0);
        Pv = TailNum / numel(DXY);
    case 'R' %(Y > X)
        TailNum = sum(DXY < 0) + 0.5*sum(DXY == 0); 
        Pv = TailNum / numel(DXY);
end

H = double(Pv <= Alpha);

if nargout >= 3
    CI = prctile(DXY, Pct);
end