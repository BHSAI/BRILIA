%labelSeqQuality will determine if a V(D)J sequence is fully annotated, is
%a productive/nonproductive VDJ, or is invalid/incomplete. The label
%results are:
%'N': Sequence has stop codon in the CDR3 or out of frame error
%'M': Sequence has stop codon in V and J. Could be sequencing error or
%     pseudogene.
%'I': Invalid or incomplete annotation, mostly likely caused by error, 
%     partial sequence, or short non-VDJ seq that resembles a junction,
%     excessive mismatches with the predicted germline sequence, or > 15
%     deletions on a V/D/J segment, or CDR length < 4
%'Y': Fully translatable to AA and all annotation values exists
%
%  VDJdata = labelSeqQuality(VDJdata, Map)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    Map: structure of locations of VDJheader
%
%  OUTPUT
%    VDJdata: VDJdata with the "Function" column filled in with N, Y, or I
%
function VDJdata = labelSeqQuality(VDJdata, Map, ThreshHold)
if isempty(VDJdata); return; end
Map = getVDJmapper(Map);
Chain = lower(Map.Chain);
if nargin < 3
    ThreshHold = 0.4;
end

for c = 1:length(Map.Chain)
    SeqIdx      = Map.([Chain(c) 'Seq']);
    RefSeqIdx   = Map.([Chain(c) 'RefSeq']);
    LengthIdx   = Map.([Chain(c) 'Length']);
    GeneNumIdx  = Map.([Chain(c) 'GeneNum']); 
    GeneNameIdx = Map.([Chain(c) 'GeneName']);
    CDR3Idx     = Map.([Chain(c) 'CDR3']);
    FunctIdx    = Map.([Chain(c) 'Funct']);
    DelIdx      = Map.([Chain(c) 'Del']);
      
    for j = 1:size(VDJdata,1)
        
        %Extract necessary information
        Seq    = VDJdata{j, SeqIdx};
        VMDNJ  = cell2mat(VDJdata(j,LengthIdx));
        VNum   = VDJdata{j, GeneNumIdx(1)};      
        VName  = VDJdata{j, GeneNameIdx(1)};
        CDR3S  = VDJdata{j, CDR3Idx(3)};
        CDR3E  = VDJdata{j, CDR3Idx(4)};
        CDR3Len = VDJdata{j, CDR3Idx(2)}; %CDR3 length is >= 5 AA and <= 30 AA including 104C and 118W
        Del    = cell2mat(VDJdata(j, DelIdx));
 
        %Make sure all necessary information is available
        VDJdata{j, FunctIdx} = 'I'; %Start with an incomplete, and then determine if it's a Y
        
        if isempty(Seq); continue; end
        if isempty(VNum); continue; end
        if isempty(VName); continue; end
        if isempty(CDR3S); continue; end
        if isempty(CDR3E); continue; end
        if isempty(CDR3Len); continue; end
        if isempty(Del); continue; end
        
        %Make sure info makes sense (I for invalid)
        %NOTE! CDR3Len must add 2 because this is for the PRE_IMGT_CDR3Length!
        if min(VMDNJ) < 0 || sum(VMDNJ) ~= length(Seq) || CDR3S < 1 || CDR3E > length(Seq) || CDR3Len < (5+2) || CDR3Len > (30+2) || max(Del) > 15
            VDJdata{j, FunctIdx} = 'I';
            continue
        end

        %Sanity check - make sure there there isn't too many mismatch 20 nts nearby in FR3 and 4, and overall alignment
        RefSeq = VDJdata{j, RefSeqIdx};
        if ~isempty(RefSeq)
            MatchLoc = cmprSeqMEX(Seq, RefSeq);
            FR3S = max(CDR3S - 20 + 1, 1);
            FR4E = min(CDR3E + 20, length(RefSeq));
            AllScore = sum(MatchLoc)/length(RefSeq);
            FR3Score = sum(MatchLoc(FR3S:CDR3S-1))/(CDR3S-FR3S);
            FR4Score = sum(MatchLoc(CDR3E+1:FR4E))/(FR4E-CDR3E);
            if min([AllScore FR3Score FR4Score]) < (1 - ThreshHold)
                VDJdata{j, FunctIdx} = 'I';
                continue
            end
        end
        
        %See if the CDR3 nt length is a multiple of 3 for in-frame junction
        if mod(CDR3E - CDR3S + 1, 3) > 0 
            VDJdata{j,FunctIdx} = 'N';
            continue
        end

        %See if there is any stop codon in the CDR3 sequence.
        AAseq = nt2aa(Seq(CDR3S:CDR3E), 'ACGTonly', false, 'frame', 1);
        if ~isempty(find(AAseq == '*', 1))
            VDJdata{j, FunctIdx} = 'N';
            continue
        end
        
        %Stop codons in V/J are classified as "M" for maybe sequencing error or not.
        AAseq = nt2aa(Seq, 'ACGTonly', false, 'frame', mod(CDR3S-1, 3) + 1);
        if ~isempty(find(AAseq == '*', 1))
            VDJdata{j, FunctIdx} = 'M';
            continue
        end
        
        %If all passes, then assume Y
        VDJdata{j, FunctIdx} = 'Y';
    end
end