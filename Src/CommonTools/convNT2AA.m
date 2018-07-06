%convNT2AA is a wrapper for nt2aa that converts nt to aa, but with extra
%checks to ensure that only valid characters are passed to nt2aa matlab
%built-in-function. Note that nt2aa('gax', 'acgtonly', false) will return
%an error, but not convNT2AA('gax', 'acgtonly', false).

function AA = convNT2AA(NT, varargin)
AmbigIdx = regexpi(NT, '[^ACGTU]');
NT(AmbigIdx) = 'N';
AA = nt2aa(NT, varargin{:});