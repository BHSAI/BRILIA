%labelNonprodVDJ will determine if a V(D)J sequence can be translated to AA
%without error. Will label the "Functional" column of VDJdata with:
%'N': Sequence has pseudogene, stop codon *, or out of frame error
%'Y': Entire sequence is translatable to AA
%'M': CDR3 sequence end or start is cutoff to determine for certain
%'I': Invalid annotation, mostly likely caused by error
%
%  VDJdata = labelNonprodVDJ(VDJdata,VDJheader,DB)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    VDJheader: main BRILIA header cell
%    DB: gene database structure(getCurrentDatabase.m)
%
%  OUTPUT
%    VDJdata: modified VDJdata with the "Function" column filled in.

function VDJdata = labelNonprodVDJ(VDJdata,VDJheader,DB)
%Determine chain
M = getMapHeaderVar(DB.MapHeader);
[H, L, Chain] = getAllHeaderVar(VDJheader);

%For each chain label N, Y, M for functionality
for k = 1:length(Chain)
    %Determine chain header locator
    if Chain(k) == 'H'
        B = H;
    else
        B = L;
    end
    
    for j = 1:size(VDJdata,1)
        %Extract necessary information
        Seq = VDJdata{j,B.SeqLoc};
        Vnum = VDJdata{j,B.GeneNumLoc(1)};      
        Vname = VDJdata{j,B.GeneNameLoc(1)};
        CDR3s = VDJdata{j,B.CDR3Loc(3)};
        CDR3e = VDJdata{j,B.CDR3Loc(4)};
        CDR3len = VDJdata{j,B.CDR3Loc(2)};
        VMDNJ = cell2mat(VDJdata(j,B.LengthLoc));
 
        %Make sure all necessary information is available
        if isempty(Seq); continue; end
        if isempty(Vnum); continue; end
        if isempty(Vname); continue; end
        if isempty(CDR3s); continue; end
        if isempty(CDR3e); continue; end
        if isempty(CDR3len); continue; end
        
        %Make sure info makes sense (I for invalid)
        if min(VMDNJ) < 0 || sum(VMDNJ) ~= length(Seq); 
            VDJdata{j,B.FunctLoc} = 'I';
            continue; 
        end
        if CDR3s < 1 || CDR3e > length(Seq); 
            VDJdata{j,B.FunctLoc} = 'I';
            continue; 
        end

        %See if the CDR3 nt length is a multiple of 3 for in-frame junction
        if mod(CDR3e - CDR3s + 1,3) > 0
            VDJdata{j,B.FunctLoc} = 'N';
            continue;
        end

        %See if CDR3 length is >= 5 AA and <= 30 AA including 104C and 118W
        if isempty(CDR3len) %Forgot to calculate?
            CDR3len = (CDR3e - CDR3s + 1)/3;
        end
        if CDR3len < 5 || CDR3len > 30
            VDJdata{j,B.FunctLoc} = 'N';
            continue
        end

        %See if there is any stop codon in the CDR3 sequence
        %ReadFrame = mod(CDR3s+2,3)+1; No longer needed, as we only check
        %from CDR3
        BadCharIdx = regexpi(Seq,'[^ACGTU]');
        Seq(BadCharIdx) = 'N';
        AAseq = nt2aa(Seq(CDR3s:CDR3e),'ACGTonly','false','frame',1);
        if ~isempty(regexp(AAseq,'\*','once'))
            VDJdata{j,B.FunctLoc} = 'N';
            continue
        end

        %See if this is a V pseudogenes
        HavePseudo = 0; %Assume maybe for now, until you confirm a "F".
        for v = 1:length(Vnum)
            if Chain(k) == 'H'
                if ~isempty(regexpi(DB.Vmap{Vnum(v),M.FunctLoc},'P','once'));
                    HavePseudo = 1;
                    break
                end
            else
                if ~isempty(regexpi(Vname,'IGKV','once')) %Kappa
                    if ~isempty(regexpi(DB.Vkmap{Vnum(v),M.FunctLoc},'P','once'));
                        HavePseudo = 1;
                        break
                    end
                else %Lambda
                    if ~isempty(regexpi(DB.Vlmap{Vnum(v),M.FunctLoc},'P','once'));
                        HavePseudo = 1;
                        break
                    end
                end
            end
        end
        if HavePseudo == 1 && length(Vnum) == 1 %Have only pseudo, so not functional
            VDJdata{j,B.FunctLoc} = 'N';
        elseif HavePseudo == 1 && length(Vnum) > 1 %Have functional and pseudo, so maybe. Should use fixDegenVDJ to remove pseudo amongst functionals.
            VDJdata{j,B.FunctLoc} = 'M';
        else %If it passes all test, it's functional
            VDJdata{j,B.FunctLoc} = 'Y';
        end
    end
end
