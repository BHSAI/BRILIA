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

function VDJdata = labelNonprodVDJ(VDJdata,Map,DB)
%Determine chain
M = getMapHeaderVar(DB.MapHeader);
%For each chain label N, Y, M for functionality
for k = 1:length(Map.Chain)
    %Determine chain header locator
    if Map.Chain(k) == 'H'
        GeneNumLoc  = Map.hGeneNum; 
        SeqLoc  = Map.hSeq;
        LengthLoc  = Map.hLength;
        GeneNameLoc  = Map.hGeneName;
        FunctLoc = Map.hFunct;
        CDR3Loc = Map.hCDR3;
    else
        GeneNumLoc  = Map.lGeneNum; 
        SeqLoc  = Map.lSeq;
        LengthLoc  = Map.lLength;
        GeneNameLoc  = Map.lGeneName;
        FunctLoc = Map.lFunct;
        CDR3Loc = Map.lCDR3;
    end
    
    for j = 1:size(VDJdata,1)
        %Extract necessary information
        Seq = VDJdata{j,SeqLoc};
        Vnum = VDJdata{j,GeneNumLoc(1)};      
        Vname = VDJdata{j,GeneNameLoc(1)};
        CDR3s = VDJdata{j,CDR3Loc(3)};
        CDR3e = VDJdata{j,CDR3Loc(4)};
        CDR3len = VDJdata{j,CDR3Loc(2)};
        VMDNJ = cell2mat(VDJdata(j,LengthLoc));
 
        %Make sure all necessary information is available
        if isempty(Seq); continue; end
        if isempty(Vnum); continue; end
        if isempty(Vname); continue; end
        if isempty(CDR3s); continue; end
        if isempty(CDR3e); continue; end
        if isempty(CDR3len); continue; end
        
        %Make sure info makes sense (I for invalid)
        if min(VMDNJ) < 0 || sum(VMDNJ) ~= length(Seq)
            VDJdata{j,FunctLoc} = 'I';
            continue
        end
        if CDR3s < 1 || CDR3e > length(Seq)
            VDJdata{j,FunctLoc} = 'I';
            continue
        end

        %See if the CDR3 nt length is a multiple of 3 for in-frame junction
        if mod(CDR3e - CDR3s + 1,3) > 0
            VDJdata{j,FunctLoc} = 'N';
            continue
        end

        %See if CDR3 length is >= 5 AA and <= 30 AA including 104C and 118W
        if isempty(CDR3len) %Forgot to calculate?
            CDR3len = (CDR3e - CDR3s + 1)/3;
        end
        if CDR3len < 5 || CDR3len > 30
            VDJdata{j,FunctLoc} = 'N';
            continue
        end

        %See if there is any stop codon in the CDR3 sequence
        %ReadFrame = mod(CDR3s+2,3)+1; No longer needed, as we only check
        %from CDR3
        BadCharIdx = regexpi(Seq,'[^ACGTU]');
        Seq(BadCharIdx) = 'N';
        AAseq = convNT2AA(Seq(CDR3s:CDR3e),'ACGTonly','false','frame',1);
        if ~isempty(regexp(AAseq,'\*','once'))
            VDJdata{j,FunctLoc} = 'N';
            continue
        end

        %See if this is a V pseudogenes
        HavePseudo = 0; %Assume maybe for now, until you confirm a "F".
        for v = 1:length(Vnum)
            if Map.Chain(k) == 'H'
                if ~isempty(regexpi(DB.Vmap{Vnum(v),M.FunctLoc},'P','once'))
                    HavePseudo = 1;
                    break
                end
            else
                if ~isempty(regexpi(Vname,'IGKV','once')) %Kappa
                    if ~isempty(regexpi(DB.Vkmap{Vnum(v),M.FunctLoc},'P','once'))
                        HavePseudo = 1;
                        break
                    end
                else %Lambda
                    if ~isempty(regexpi(DB.Vlmap{Vnum(v),M.FunctLoc},'P','once'))
                        HavePseudo = 1;
                        break
                    end
                end
            end
        end
        if HavePseudo == 1 && length(Vnum) == 1 %Have only pseudo, so not functional
            VDJdata{j,FunctLoc} = 'N';
        elseif HavePseudo == 1 && length(Vnum) > 1 %Have functional and pseudo, so maybe. Should use fixDegenVDJ to remove pseudo amongst functionals.
            VDJdata{j,FunctLoc} = 'M';
        else %If it passes all test, it's functional
            VDJdata{j,FunctLoc} = 'Y';
        end
    end
end
