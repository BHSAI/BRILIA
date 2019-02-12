%calcBhattStats will calculate the Bhattacharya distance and coefficient
%for comparing two discrete distibutions.
%
%  [Dist, Coef] = calcBhattStat(Y1, Y2)
%
%  INPUTS
%    Y1: Mx1 matrix of frequencies of dataset 1
%    Y2: Mx1 matrix of frequencies of dataset 2
%
%  OUTPUTS
%    Dist: Bhattacharya distance
%    Coef: Bhattacharya coefficient

function [Dist, Coef] = calcBhattStat(Y1, Y2)
Y1 = Y1/sum(Y1);
Y2 = Y2/sum(Y2);
Coef = sum((Y1.*Y2).^(0.5));
Dist = -log(Coef);