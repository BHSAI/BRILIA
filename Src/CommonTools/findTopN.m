%findTopN will find the values and index of the top N values in a vector.
%
%  Y = findTopN(X, N)
%
%  [Y, Idx] = findTopN(X, N)
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
%    [Y, Idx] = findTopN(X, 3)
%    Y =
%         10     9     8
%    Idx =
%         10     9     8


function [Y, Idx] = findTopN(X, N)
[Y, Idx] = sort(X, 'descend');
if N > length(X) 
    return
else
    Y = Y(1:N);
    Idx = Idx(1:N);    
end