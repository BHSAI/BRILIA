%binormPDF will switch between using either a binomial or a normal
%distribution to return the true or estimated probability of X successes
%out of N trials given a success probability P. 
%
%  Y = binormPDF(X, N, P)
%
%  Y = binormPDF(X, N, P, Method)
%
%  INPUT
%    X: number of successes
%    N: total number of trials
%    P: probability of success
%    Method ['b' 'n' 'a']: force use of binomial(b), normal(n), or any(a).
%       With "any", normal approx. is used if N > 50 and N*P*(1-P) > 10.
%       The binomial equation is used otherwise.
%
%  OUTPUT
%    Y: binomial or approximated normal probability density
%
%  EXAMPLE
%    X = 1:50;
%    N = 100;
%    P = 0.5;
%    Y1 = binormPDF(X, N, P, 'bi');
%    Y2 = binormPDF(X, N, P, 'norm');

function Y = binormPDF(X, N, P, varargin)
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
if Method == 'a' && N > 50 && N*P*(1-P) > 10
    UseNorm = 1;
else
    UseNorm = 0;
end

if Method == 'b' || (Method == 'a' && ~UseNorm)
    Y = binopdf(X, N, P);
elseif Method == 'n' || (Method == 'a' && UseNorm)
    Mu = N*P;
    Sig = sqrt(N*P*(1-P));
    Y = normcdf(X+0.5, Mu, Sig) - normcdf(X-0.5, Mu, Sig);
end