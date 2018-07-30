%findCDR1 will return the CDR1 information, such as AA, AA length, and
%start and end location of this region for the Sequence in VDJdata.
%
%  VDJdata = findCDR1(VDJdata, Map, DB)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    Map: structure of indeix BRILIA data columns
%    DB: Gene database structure (getCurrentDatabase.m)
%
%  OUTPUT
%    VDJdata: modified VDJdata with the CDR information filled
%
%  NOTE
%    This must be run after finding V(D)J annotations.
%
function VDJdata = findCDR1(VDJdata, Map, DB, varargin)
M = getMapHeaderVar(DB.MapHeader);
for c = 1:length(Map.Chain)
    C = lower(Map.Chain(c));
    DelIdx = Map.([C 'Del'])(1);
    SeqIdx = Map.([C 'Seq']);
    LengthIdx = Map.([C 'Length'])(1);
    GeneNumIdx = Map.([C 'GeneNum'])(1); 
    GeneNameIdx = Map.([C 'GeneName'])(1);
    ValidIdx = find(~any(cellfun('isempty', VDJdata(:, [DelIdx; SeqIdx; LengthIdx; GeneNumIdx; GeneNameIdx])), 2));
    CDRIdx = Map.([C 'CDR1']);
    
    for k = 1:length(ValidIdx)
        j = ValidIdx(k);
        %Get basic annotation information
        Seq = VDJdata{j, SeqIdx};
        Vlen = VDJdata{j, LengthIdx(1)};
        Vdel = VDJdata{j, DelIdx(1)};
        VmapNum = VDJdata{j, GeneNumIdx(1)}(1);

        %Determine the Vxmap and CDRX location
        if strcmpi(C, 'H') %Heavy chain
            MapName = 'Vmap';
        else %Light chain, determine locus Vxmap
            if contains(VDJdata{j, GeneNameIdx}, 'IGKV')
                MapName = 'Vkmap';
            else
                MapName = 'Vlmap';
            end
        end
        
        %Get basic ref gene information
        if any(cellfun('isempty', DB.(MapName)(VmapNum, [M.Seq M.CDR1s M.CDR1e]))); continue; end
        CDR1SeqLen = length(DB.(MapName){VmapNum, M.Seq});
        CDR1sRef = DB.(MapName){VmapNum, M.CDR1s};
        CDR1eRef = DB.(MapName){VmapNum, M.CDR1e};
        if min([CDR1sRef CDR1eRef]) == 0; continue; end
        CDR1RefSeq = nt2aa(DB.(MapName){VmapNum, M.Seq}(CDR1sRef:CDR1eRef), 'ACGTonly', false);
        
        %Compute location of CDRX from Seq
        Vshift = Vlen + Vdel - CDR1SeqLen;
        CDR1s = Vshift + CDR1sRef;
        CDR1e = Vshift + CDR1eRef;
        
        if min([CDR1s CDR1e]) > 0 
            CDR1Seq = nt2aa(Seq(CDR1s:CDR1e), 'ACGTonly', false);
            CDR1Seq(CDR1Seq ~= CDR1RefSeq) = lower(CDR1Seq(CDR1Seq ~= CDR1RefSeq));
        else
            CDR1Seq = CDR1RefSeq;
        end

        VDJdata(j, CDRIdx) = {CDR1Seq length(CDR1Seq) CDR1s CDR1e};
    end
end