%convNT2AA is a wrapper for nt2aa for preventing a sequence from being
%treated as a start codon when it's midway. For instance, 'TTG' is either
%an 'M' if it's starting a sequence, or an 'L' if it's in the middle.
%
%  INPUT
%    NT: cell array of char or char array of DNA/RNA sequence(s)
%    varargin: additional inputs for the regular nt2aa function
%
%  OUTPUT
%    AA: cell array of char or char array of amino acid sequence(s)
%
%  EXAMPLE
%    NT = 'TTGTTG'
%    %Regular nt2aa will treat start TTG as M, not L
%    AA_nt2aa = nt2aa(NT)
%    AA_nt2aa = 
%       'ML'
%    %This function will treat TTG as non-start codon
%    AA_convNT2AA = convNT2AA(NT)
%    AA_convNT2AA = 
%       'LL'
%
%    %Works on cell input
%    NT = {NT; NT}
%    AA_convNT2AA = 
%       'LL'
%       'LL'
%

function AA = convNT2AA(NT, varargin)
error('%s: This is obsolete.', mfilename)

NT = strcat('TTG', NT);
AA = nt2aa(NT, varargin{:});
if iscell(AA)
    parfor j = 1:numel(AA)
        AA{j} = AA{j}(2:end);
    end
else
    AA = AA(2:end);
end
