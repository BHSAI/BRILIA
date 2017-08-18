%findCDR2 will return the CDR2 information, such as AA, AA length, and
%start and end location of this region for the Sequence in VDJdata.
%
%  VDJdata = findCDR2(VDJdata, VDJheader, DB)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    VDJheader: main BRILIA header cell
%    DB: Gene database structure (getCurrentDatabase.m)
%
%  OUTPUT
%    VDJdata: modified VDJdata with the CDR2 information filled
%
%  NOTE
%    This must be run after finding V(D)J annotations.
%
%  See also findCDR1, findCDR3

function VDJdata = findCDR2(VDJdata,VDJheader,DB)
%Determine chain and extract key locations
M = getMapHeaderVar(DB.MapHeader);
H = getHeavyHeaderVar(VDJheader);
L = getLightHeaderVar(VDJheader);
if H.SeqLoc > 0 && L.SeqLoc > 0
    Chain = 'HL';
elseif H.SeqLoc > 0
    Chain = 'H';
elseif L.SeqLoc > 0
    Chain = 'L';
end

for k = 1:length(Chain)
    %Decide what header locator to use
    if Chain(k) == 'H';
        B = H;
    else
        B = L;
    end
    
    for j = 1:size(VDJdata,1)
        %Get basic annotation information
        Seq = strrep(VDJdata{j,B.SeqLoc},'X','N');
        Vlen = VDJdata{j,B.LengthLoc(1)};
        Vdel = VDJdata{j,B.DelLoc(1)};
        VmapNum = VDJdata{j,B.GeneNumLoc(1)};

        %Make sure all info is provided for CDRX calculation 
        if isempty(Seq); continue; end
        if isempty(Vlen); continue; end
        if isempty(Vdel); continue; end
        if isempty(VmapNum); continue; end
        
        %Determine the Vxmap and CDRX location
        if Chain(k) == 'H' %Heavy chain
            Vxmap = DB.Vmap;
        else %Light chain, determine locus Vxmap
            %Determine the locus and anchor loc
            Vname = VDJdata{j,L.GeneNameLoc(1)};
            if ~isempty(regexpi(Vname,'IGKV','once')) %Kappa
                Vxmap = DB.Vkmap;
            else %Lambda
                Vxmap = DB.Vlmap;
            end
        end
        
        %Get basic ref gene information
        Vseq = Vxmap{VmapNum,M.SeqLoc};
        CDR2sRef = Vxmap{VmapNum,M.CDR2sLoc};
        CDR2eRef = Vxmap{VmapNum,M.CDR2eLoc};

        %Make sure all info is provided for ref CDRX calculation 
        if isempty(Vseq); continue; end
        if isempty(CDR2sRef) || CDR2sRef == 0; continue; end
        if isempty(CDR2eRef) || CDR2eRef == 0; continue; end
        
        %Compute location of CDRX from Vseq
        CDR2aaRef = upper(nt2aa(Vseq(CDR2sRef:CDR2eRef),'ACGTonly','false'));
        
        %Compute location of CDRX from Seq
        Vshift = Vlen + Vdel - length(Vseq);
        CDR2s = Vshift + CDR2sRef;
        CDR2e = Vshift + CDR2eRef;
        if min([CDR2s CDR2e]) <= 0 %Need to CDRX of the V ref seq
            CDR2aa = CDR2aaRef;
        else %Need to get CDRX of the real seq
            CDR2aa = nt2aa(Seq(CDR2s:CDR2e),'ACGTonly','false');
            CDR2aa(CDR2aa ~= CDR2aaRef) = lower(CDR2aa(CDR2aa ~= CDR2aaRef)); %Make SHMs lowercase
        end

        %Fill in the VDJdata
        VDJdata(j,B.CDR2Loc) = {CDR2aa length(CDR2aa) CDR2s CDR2e};
    end
end
