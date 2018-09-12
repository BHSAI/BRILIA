%findBotN will find the values and index of the top N values in a vector.
%
%  Y = findBotN(X, N)
%
%  [Y, Idx] = findBotN(X, N)
%
%  INPUT
%    X: vector of values
%    N: the N number of top values to collect from X
%
%  OUTPUT
%    Y: the top N number of values sorted in descending order
%    Idx: the location of Y in X   Y = X(Idx)
%
%  EXAMPLE
%    X = [1 2 3 4 5 6 7 8 9 10]
%    [Y, Idx] = findBotN(X, 3)
%    Y =
%         1     2     3
%    Idx =
%         1     2     3
%

function [Y, Idx] = findBotN(X, N)
[Y, Idx] = sort(X, 'ascend');
if N > length(X) 
    return
else
    Y = Y(1:N);
    Idx = Idx(1:N);    
end