%labelNonprodVDJ will look for sequences that have stop codons,
%out-of-frame CDR3s, pseudo gene, and then label a "N" under the
%"Functional" column of VDJdata for non-funcitonal. All others with
%conserved C and W of CDR3 regions will be labeled 'Y'.
%
%  VDJdata = labelNonprodVDJ(VDJdata,NewHeader)
%
%  VDJdata = labelNonprodVDJ(VDJdata,NewHeader,Vmap,Dmap,Jmap)

function VDJdata = labelNonprodVDJ(VDJdata,NewHeader,varargin)
%Extract the VDJ database
if length(varargin) >= 3
    [Vmap,Dmap,Jmap] = deal(varargin{:});
else
    [Vmap, ~, ~] = getCurrentDatabase;
end
getHeaderVar;

for j = 1:size(VDJdata,1)
    try
        %See if the CDR3start to end is multiple of 3 for in-frame junction
        CDR3start = VDJdata{j,CDR3Loc(3)};
        CDR3end = VDJdata{j,CDR3Loc(4)};
        if mod(CDR3end - CDR3start + 1,3) > 0
            VDJdata{j,FunctLoc} = 'N';
            continue;
        end
        
        %See if CDR3 length makes sense, must be longer than 4 AA,
        %including 104C and 118W
        CDR3Len = VDJdata{j,CDR3Loc(2)};
        if isempty(CDR3Len) %Forgot to calculate?
            CDR3Len = (CDR3end - CDR3start + 1)/3;
        end
        if CDR3Len <= 4
            VDJdata{j,FunctLoc} = 'N';
            continue
        end

        %See if there are any stop codons in the sequence
        Seq = VDJdata{j,SeqLoc};
        ReadFrame = mod(CDR3start+2,3)+1;
        if min(isempty(Seq)) == 1 || isempty(ReadFrame)
            continue
        end
        BadCharIdx = regexpi(Seq,'[^ACGTU]');
        Seq(BadCharIdx) = 'N';
        AAseq = nt2aa(Seq,'ACGTonly','false','frame',ReadFrame);
        if ~isempty(regexp(AAseq,'\*','once'))
            VDJdata{j,FunctLoc} = 'N';
            continue
        end

        %See if there are any pseudogenes
        Vnum = VDJdata{j,FamNumLoc(1)};      
        HavePseudo = 0; %Assume maybe for now, until you confirm a "F".
        for v = 1:length(Vnum)
            if ~isempty(regexpi(Vmap{Vnum(v),7},'P','once'));
                HavePseudo = 1;
                break
            end
        end
        if HavePseudo == 1 && length(Vnum) == 1 %Have only pseudo, so not functional
            VDJdata{j,FunctLoc} = 'N';
        elseif HavePseudo == 1 && length(Vnum) > 1 %Have functional and pseudo, so maybe. Should use fixDegenVDJ to remove pseudo amongst functionals.
            VDJdata{j,FunctLoc} = 'M';
        else %If it passes all test, it's functional
            VDJdata{j,FunctLoc} = 'Y';
        end
    catch
        ErrorMsg = sprintf('Errored at %s, sequence # %d',mfilename,j);
        disp(ErrorMsg);
        VDJdata{j,FunctLoc} = 'E';
        VDJdata{j,MiscLoc} = ErrorMsg;
    end
end