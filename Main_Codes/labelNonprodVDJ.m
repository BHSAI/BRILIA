%labelNonprodVDJ will look for sequences that have stop codons,
%out-of-frame CDR3s, and then label a "N" under the "Functional" column of
%VDJdata for non-funcitonal. All others with conserved C and W of CDR3
%regions will be labeled 'Y'.

function VDJdata = labelNonprodVDJ(VDJdata,NewHeader,varargin)
%Extract the VDJ database
if length(varargin) == 3
    Vmap = varargin{1};
    Dmap = varargin{2};
    Jmap = varargin{3};
else
    [Vmap, ~, ~] = getCurrentDatabase;
end
getHeaderVar;

for j = 1:size(VDJdata,1)
    %See if CDR3 exist
    CDR3Len = VDJdata{j,CDR3Loc(2)};
    if CDR3Len <= 4
        VDJdata{j,FunctLoc} = 'N';
        continue
    end
    Seq = VDJdata{j,SeqLoc};
    CDR3start = VDJdata{j,CDR3Loc(3)}+2;
    
    %See if there are any stop codons in the sequence
    ReadFrame = mod(CDR3start,3)+1;
    if min(isempty(Seq)) == 1 || isempty(ReadFrame)
        continue
    end
    AAseq = nt2aa(Seq,'ACGTonly','false','frame',ReadFrame);
    if ~isempty(regexp(AAseq,'\*','once'))
        VDJdata{j,FunctLoc} = 'N';
        continue
    end
    
    %See if there are any pseudogenes
    Vnum = VDJdata{j,FamNumLoc(1)};
    Vstatus = 0; %Assume maybe for now, until you confirm a "F".
    for v = 1:length(Vnum)
        if ~isempty(regexpi(Vmap{Vnum(v),7},'P','once'));
            Vstatus = 1;
            break
        end
    end
    if Vstatus == 1
        VDJdata{j,FunctLoc} = 'M';
        continue
    end
    
    %If it passes all test, it's functional
    VDJdata{j,FunctLoc} = 'Y';
end