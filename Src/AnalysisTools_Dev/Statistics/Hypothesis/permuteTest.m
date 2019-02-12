%permuteTest will perform a permutation test to see if two vectors of data
%have the same or different means. This is an absolute test that should
%outperform a ttest, but the issue is that the permutations can easily
%exceed a reasonable number. 
%
%  INPUT
%    X: vector of values from sample A
%    Y: vector of values from sample B
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
%    H: 1 for reject H0, 0 for fail to reject H0
%    P: p-value for observing difference in means of X and Y
%    CI: 1x2 matrix for the (1-Alpha)*100% confidence interval of the
%      difference in medians
%    DM: Difference in mean of X and y
%
%  EXAMPLE
%    X = [2 3 4 5 3 2];
%    Y = [8 7 10 11 7];
%    Alpha = 0.05;
%    [H, P, CI] = permuteTest(X, Y, 'Alpha', Alpha)
%    H =
%       1
%    P =
%       0.0023
%    CI =
%      -3.5667    3.5667
% 
%    [H, P, CI] = ttest2(X, Y, 'alpha', Alpha)
%    H =
%        1
%    P =
%        1.9901e-04
%    CI =
%       -7.4770   -3.3896
% 
function [H, P, CI, DM] = permuteTest(X, Y, varargin)
P = inputParser;
addParameter(P, 'Alpha', 0.05, @(x) isnumeric(x) && x >= 0 && x <= 1);
addParameter(P, 'Tail', 'B', @(x) startsWith(upper(x), {'L', 'R', 'B'}));
addParameter(P, 'Iter', 100000, @(x) x >= 0 && mod(x, 1) == 0); %factorial(9) is max, so go no less than that
parse(P, varargin{:});
Alpha = P.Results.Alpha;
Tail = P.Results.Tail;
Iter = P.Results.Iter;

DMobs = mean(Y) - mean(X);
if (numel(X) + numel(Y)) > 8
    warning('%s: Did not code for N > 8 permutation tests. Using bootstrap instead.', mfilename);
    XY = [X(:); Y(:)];
    PermIdx = randi(numel(XY), Iter, numel(XY));
    AllPerm = XY(PermIdx);
else
    AllPerm = perms([X(:); Y(:)]);
end

Mx = mean(AllPerm(:, 1:numel(X)), 2);
My = mean(AllPerm(:, numel(X)+1:end), 2);
DM = My-Mx;

switch upper(Tail)
    case 'B' %2-tail both end
        Pct = [Alpha/2 1-Alpha/2]*100;
        P = sum(DM >= DMobs | DM <= -DMobs)/numel(DM);
    case 'L' %Left tail
        Pct = [Alpha 1]*100;
        P = sum(DM >= DMobs)/numel(DM);
    case 'R' %Right tail
        Pct = [0 1-Alpha]*100;
        P = sum(DM <= -DMobs)/numel(DM);
end
CI = prctile(DM, Pct);

H = double(P <= Alpha);