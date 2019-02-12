function [ChiSq, PValue] = calcChiStat(Y1, Y2, Option)
Loc = Y1 ~= 0 & Y2 ~= 0;
if nargin < 3 || startsWith(Option, 'norm', 'IgnoreCase', true)
    R = Y1(Loc)/sum(Y1(Loc));
    S = Y2(Loc)/sum(Y2(Loc));
else
    R = Y1(Loc);
    S = Y2(Loc);
end

SumR = sum(R);
SumS = sum(S);
ChiSq = sum( (sqrt(SumS/SumR)*R - sqrt(SumR/SumS)*S).^2 ./ (R + S) );

PValue = chi2cdf(ChiSq, 1);