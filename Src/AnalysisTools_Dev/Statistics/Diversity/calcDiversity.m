%calcDiversity will compute the true diversity given an entropy name.
%
%  [D, I] = calcDiversity(Option, X, Q)
%
%  INPUT
%    Option: the entropy name
%      richness
%      shannon
%      simpson
%      ginisimpson
%      hcdt
%      renyi
%    X: frequency data
%    Q: diversity order
%
%  OUTPUT
%    D: true diversity
%    I: entropy index

function [D, I] = calcDiversity(Option, X, Q)
if nargin < 3
     Q = 1;
end

switch lower(Option)
    case 'richness'
        I = numel(X);
        D = I;
    case 'shannon'
        P = X(:)/sum(X(:));
        I = -sum(P.*log(P));
        D = exp(I);
    case 'simpson'
        P = X(:)/sum(X(:));
        I = sum(P.^2);
        D = 1/I;
    case {'ginisimpson', 'gini-simpson'}
        P = X(:)/sum(X(:));
        I = 1 - sum(P.^2);
        D = 1/(1-I);
    case 'hcdt'
        P = X(:)/sum(X(:));
        I = (1 - sum(P.^Q)) / (Q - 1);
        D = ((1 - (Q - 1)*I))^(1/(1-Q));
    case 'renyi'
        P = X(:)/sum(X(:));
        I = (-log(sum(P.^Q))) / (Q - 1);
        D = exp(I);
end