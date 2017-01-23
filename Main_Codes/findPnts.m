%findPnts will take a sequence and find p-nucleotides.
%
%  Ploc = findPnts(RefSeq, NonPrange, Pside)
%
%  Ploc = findPnts(RefSeq, NonPrange, Pside, RefDelLeft, RefDelRight)
%
%  INPUT
%    RefSeq: nucleotide sequence of the germline sequence
%    NonPrange [s1,s2]: index range of non-Pnucleotides
%    Pside ['left','right','both']: find P nts left or right of range values
%    RefDelLeft: Integer >= 0, will only find P if RefDelLeft = 0;
%    RefDelRight: Integer >= 0, will only find P if RefDelRight = 0;
%
%  NOTE
%    RefDelLeft and RefDelRight must be 0 for P-nts to be found. The reason
%    for including these inputs is in case you want to use findPnts in a
%    function, but do not want to write extra code to check if P nts can
%    exist.
%
%  EXAMPLE
%    RefSeq = 'ACGTTATGCATGAA'
%    NonPrange = [1,8];
%    Pside = 'right';
%    Ploc = findPnts(RefSeq,NonPrange,Pside)
%    Ploc = 
%      0  0  0  0  0  0  0  0  1  1  1  0  0  0

function Ploc = findPnts(RefSeq,NonPrange,Pside,varargin)
%Check for variable inputs
RefDelLeft = 0;
RefDelRight = 0;
if length(varargin) >= 1
    RefDelLeft = varargin{1};
    if length(varargin) >= 2
        RefDelRight = varargin{2};
    end
end

Ploc = zeros(1,length(RefSeq),'logical');

%Find P-nt left of RefSeq's NonPrange.
if RefDelLeft == 0
    if strcmpi(Pside,'left') || strcmpi(Pside,'both')
        StartLoc = NonPrange(1);
        EvalLoc = StartLoc - 1; 
        while EvalLoc >= 1
            if strcmpi(seqcomplement(RefSeq(StartLoc)),RefSeq(EvalLoc))
                Ploc(EvalLoc) = 1;
                EvalLoc = EvalLoc - 1;
                StartLoc = StartLoc + 1;
                if StartLoc > length(RefSeq)
                    break
                end
            else
                break
            end
        end
    end
end

%Find P-nt right of RefSeq's NonPrange.
if RefDelRight == 0
    if strcmpi(Pside,'right') || strcmpi(Pside,'both')
        StartLoc = NonPrange(end);
        EvalLoc = StartLoc + 1; 
        while EvalLoc <= length(RefSeq)
            if strcmpi(seqcomplement(RefSeq(StartLoc)),RefSeq(EvalLoc))
                Ploc(EvalLoc) = 1;
                EvalLoc = EvalLoc + 1;
                StartLoc = StartLoc - 1;
                if StartLoc < 1
                    break
                end
            else
                break
            end
        end    
    end
end