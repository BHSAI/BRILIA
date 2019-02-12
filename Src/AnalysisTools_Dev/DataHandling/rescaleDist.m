%rescaleDist will scale a discrete distribution such that MOST of the bars
%are similar to each other, as determined by a linear scaling factor that
%tries to minimize the differences between two bar heights.
%
%  [Yscaled, S, Err] = rescaleDist(Y, Yref, Tol)
%
%  INPUT
%    Y: MxN matrix to scale each column by Yref
%    Yref: MxV matrix of frequencies to use via a mean(Yref, 2), to rescale 
%      dist Y by per each column.
%    Tol: Tolerance in scaling factor S. Default is 0.0001.
%
%  OUTPUT
%    Yscaled: the rescaled Y frequency distribution
%    S: 1xN matrix of scaling factors such that Y*S = Yscaled
%    Err: 1xN matrix of the sum(abs(Yscaled-Y))^0.5)
%
%  EXAMPLE 1
%    Y    = [0 2 4 6 2 5 2 20]';
%    Yref = [0 1 2 3 1 5 1 10]';
%    Tol  = 1E-4;
%    [Yscaled, S, Err] = rescaleDist(Y, Yref, Tol)
%    Yscaled =
%              0
%         1.0000
%         2.0000
%         3.0000
%         1.0000
%         2.5000   %< this is the different one
%         1.0000
%        10.0000
%    S = 0.5
%    Err = 1.5811
%
%  EXAMPLE 2
%    C = {10 20  50; 
%         20 40 100;
%         10 10  50};
%    Group = [1 1 2];
%    [ScaleDist, S, Err] = rescaleDist(C, C(:, Group == 1))
%    ScaleDist =
%         [15]    [    15]    [15]
%         [30]    [    30]    [30]
%         [15]    [7.5000]    [15]
%    S =
%        1.5000    0.7500    0.3000
%    Err = 
%        2.2361    1.5811    2.2361
%
function [Y, Scalor, Err] = rescaleDist(Y, Yref, Tol)
if nargin < 3
    Tol = 1E-4;
end

IsCell = iscell(Y); %Used later to convert output to cell or not
if IsCell
    Y = cell2mat(Y);
end
if iscell(Yref)
    Yref = cell2mat(Yref);
end
if size(Y, 1) ~= size(Yref, 1)
    error('%s: Y and Yref must have same number of rows', mfilename);
end
if size(Yref, 2) > 1
    Yref = mean(Yref, 2);
end

Err = zeros(1, size(Y, 2));
Scalor = zeros(1, size(Y, 2));
for c = 1:size(Y, 2)
    [Y(:, c), Scalor(c), Err(c)] = rescaleDistPerCol(Y(:, c), Yref, Tol);
end

if IsCell
    Y = num2cell(Y);
end

%Find the rescaled distribution
function [Y, Scalor, Err] = rescaleDistPerCol(Y, Yref, Tol)
SRange = Yref./Y;
SRange = SRange(~(isnan(SRange) | isinf(SRange)));
if isempty(SRange)
    Scalor = 1;
    Err = 0;
    return
end

[Scalor, LowS, HighS, Err] = findScalor(Y, Yref, min(SRange), max(SRange));
CurErr = Inf;
while abs(CurErr - Err) > Tol
    CurErr = Err;
    [Scalor, LowS, HighS, Err] = findScalor(Y, Yref, LowS, HighS);    
end
Y = Y*Scalor;

%find the multiplication factor
function [S, LowS, HighS, MinErr] = findScalor(Y, Yref, LowS, HighS)
N = linspace(LowS, HighS, 1000);
Err = sum(abs(Y*N - Yref).^0.5);
MinErr = min(Err);
MinIdx = find(Err == MinErr, 1);
S = N(MinIdx);

if MinIdx > 1
    LowS = N(MinIdx - 1);
else
    LowS = N(1);
end

if MinIdx < length(N)
    HighS = N(MinIdx + 1);
else
    HighS = N(end);
end