%winsorize will winsorize a vector. Vector will be unsorted and P
%percentiles of both edes will be converted to the edge values.
%
%  W = winsorize(X)
%
%  W = winsorize(X, P)
%
%  INPUT
%    X: 1xN or Mx1 vector of values
%    P: scalar or 1x2 matrix percentage to convert off the edges 
%       Default is 20% from left, 20% from right)
%       Sum of P must be < 100
%
%  OUTPUT
%    W: winsorized and sorted vector of X
%    WinLoc: winsorized number logical index
%
%  EXAMPLE
%    X = 1:10
%    [W, WinLoc] = winsorize(X, 20)
%    W = 
%       3   3   3   4   5   6   7   8   8   8
%    WinLoc = 
%       1   1   0   0   0   0   0   0   1   1

function [X, WinLoc] = winsorize(X, P)

if ~isvector(X)
   error('%s: X must be a vector', mfilename);
end  
if nargin < 2
   P = [20 80];
end
if nargin == 2 
    if length(P) > 3
        error('%s: P must be a value or 1x2 matrix of percentile values', mfilename);
    elseif length(P) == 1
        P = [P 100-P];
    end
    if min(P) < 0 || max(P) >  100
        error('%s: P must have values between 0 and 100', mfilename);
    end
    if sum(P) > 100
    	error('%s: The sum of P cannot exceed 100', mfilename);
    end  
end

P  = prctile(X, P);
LocL = X < P(1);
ValL = min(X(~LocL));
LocH = X > P(2);
ValH = max(X(~LocH));
X(LocL) = ValL;
X(LocH) = ValH;

if nargout >= 2
    WinLoc = LocL | LocH;
end