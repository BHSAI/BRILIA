%calcCohenD will calculate the Cohen's D distance between two mean for size
%effect measurements. 
%
%  D = calcCohenD(X1, X2)
%
%  INPUT
%    X1: MxN matrix, where each col is a measurement, and each row is a
%        measurement type. 
%    X2: Same as X1, but for a different group of samples
%
%  OUTPUT
%    D: Mx1 matrix of absolute Cohen distance between two means for each
%      Mth measure
%
%  NOTE
%    WARNING: ALWAYS use inputs where where col = trial # and 
%      row = measurement of certain property.
%
%  EXAMPLE
%    X1 = [ 3 4 3 2; 2 3 1 2; -1  -3  -4  -5];
%    X2 = [10 7 6 9; 1 1 7 3; -4 -10 -15 -12];
%    D = calcCohenD(X1, X2)
%    D =
%       3.5355
%       0.4804
%       2.0000
%
function D = calcCohenD(X1, X2)
if iscell(X1)
    X1 = cell2mat(X1);
end
if iscell(X2)
    X2 = cell2mat(X2);
end
N1 = size(X1, 2);
N2 = size(X2, 2);
M1 = mean(X1, 2);
M2 = mean(X2, 2);
V1 = var(X1, 0, 2); %Unbiased estimator for variance w/ Bessel's correction
V2 = var(X2, 0, 2); 
STD = sqrt( ((N1-1)*V1 + (N2-1)*V2) / (N1 + N2 - 2) );
D = abs((M1-M2)./STD);