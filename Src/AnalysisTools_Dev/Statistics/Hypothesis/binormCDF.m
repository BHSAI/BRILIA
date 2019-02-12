%binormCDF will use either a binomial or a normal distribution to estimate
%the X success of N trials given a P probability. The normal distribution
%estimating is used when N >= 50.
%
%  Y = binormCDF(X, N, P)
%
%  Y = binormCDF(X, N, P, Method)
%
%  Y = binormCDF(X, N, P, ..., 'upper')
%
%  INPUT
%    X: number of successes
%    N: total number of trials
%    P: probability of success
%    Method ['b' 'n' 'a']: force use of binomial(b), normal(b), or any(a).
%       With "any", normal approx. is used if N > 50 and N*P*(1-P) > 10,
%       and binomial is used otherwise.
%    'upper': takes the upper tail of distribution, past the Xth success.
%             To get the current success and higher, use 
%             binormCDF(X-1, N, P, 'upper')
%
%  OUTPUT
%    Y: binomial or approximated normal cumulative probability function
%
%  EXAMPLE
%    X = 1:50;
%    N = 100;
%    P = 0.5;
%    Y1 = binormCDF(X, N, P, 'bi', 'upper');
%    Y2 = binormCDF(X, N, P, 'norm', 'upper');

function Y = binormCDF(X, N, P, varargin)
UpperLoc = contains(varargin, 'upper', 'ignorecase', true);
if any(UpperLoc)
    UpperCell = {'upper'};
    varargin = varargin(~UpperLoc);
else
    UpperCell = {};
end

if ~isempty(varargin)
    switch lower(varargin{1}(1))
        case 'b' %binomial
            Method = 'b';
        case 'n' %normal
            Method = 'n';
        case 'a' %any
            Method = 'a';
        otherwise
            error('%s: Unrecognized method "%s".', mfilename, varargin{1});
    end
else
    Method = 'a';
end

%Decide to use normal approx vs binomial
if Method == 'a' && N > 50 && N*P*(1-P) >= 10
    UseNorm = 1;
else
    UseNorm = 0;
end

if Method == 'b' || (Method == 'a' && ~UseNorm)
    Y = binocdf(X, N, P, UpperCell{:});
elseif Method == 'n' || (Method == 'a' && UseNorm)
    Mu = N*P;
    Sig = sqrt(N*P*(1-P));
    Y = normcdf(X+0.5, Mu, Sig, UpperCell{:});
end