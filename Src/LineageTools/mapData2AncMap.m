%mapData2AncMap will take Tdata (subset of VDJdata for same clonotype) and
%it corresponding AncMap, and then fill in RefSeq as the parent sequence of
%each child sequence.
%
%  Tdata = mapData2AncMap(Tdata, VDJheader, AncMap)
%
%  INPUT
%    Tdata: main BRILIA data cell, but selecting for only sequences that
%      belong in the same crude cluster (see clusterGene.m) and have the
%      same sequence lengths (padtrimSeqGroup.m).
%    VDJheader: main BRILIA header cell
%    AncMap: Ancestrap map matrix (calcAncMap.m or calcRootedAncMap.m)
%
%  OUTPUT
%    Tdata: modified Tdata where the RefSeq is the parent sequence of each
%      child sequence. 
%
%  NOTE
%    The RefSeq in INPUT Tdata for the 1st row should be germline RefSeqs.
%    
%    The RefSeq of the 1st row of Tdata is not touched, along with any
%      AncMap(j, 2) = 0, which is the germline parent 0.
function Tdata = mapData2AncMap(Tdata, VDJheader, AncMap)
%Determine chain and mutation locations
[H, L, Chain] = getAllHeaderVar(VDJheader);
if strcmpi(Chain, 'HL')
    SeqLoc = [H.SeqLoc; L.SeqLoc];
    RefSeqLoc = [H.RefSeqLoc; L.RefSeqLoc];
elseif strcmpi(Chain, 'H')
    SeqLoc = H.SeqLoc;
    RefSeqLoc = H.RefSeqLoc;
elseif strcmpi(Chain, 'L')
    SeqLoc = L.SeqLoc;
    RefSeqLoc = L.RefSeqLoc;
end

%Make sure AncMap is relative to Tdata position. 
AncMap = renumberAncMap(AncMap);

%For each child seq, you want the RefSeq to be the Seq of the parent.
for j = 2:size(Tdata, 1)
    if AncMap(j, 2) == 0; continue; end %Always leave root alone
    for c = 1:length(RefSeqLoc)
        Tdata(AncMap(j, 1), RefSeqLoc(c)) = Tdata(AncMap(j, 2), SeqLoc(c));
    end
end
