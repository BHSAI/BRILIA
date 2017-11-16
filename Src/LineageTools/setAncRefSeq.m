%setAncRefSeq will take Tdata (subset of VDJdata for one clonotype) and
%its corresponding AncMap, and fill in RefSeq as the parent sequence of
%each child sequence.
%
%  Tdata = setAncRefSeq(Tdata, VDJheader, AncMap)
%
%  INPUT
%    Tdata: subset of VDJdata containing one linked cluster (see
%      clusterGene.m) and have the same sequence lengths
%      (padtrimSeqGroup.m).
%    VDJheader: main BRILIA header cell
%    AncMap: ancestral map matrix (calcAncMap.m or calcRootedAncMap.m)
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

function Tdata = setAncRefSeq(Tdata, Map, AncMap)
if size(Tdata, 1) ~= AncMap
    warning('%s: Tdata and AncMap row dimension cannot be different.', mfilename);
end

%Determine chain and mutation locations
if strcmpi(Map.Chain, 'HL')
    SeqLoc = [Map.hSeq; Map.lSeq];
    RefSeqLoc = [Map.hRefSeq; Map.lRefSeq];
elseif strcmpi(Map.Chain, 'H')
    SeqLoc = Map.hSeq;
    RefSeqLoc = Map.hRefSeq;
elseif strcmpi(Map.Chain, 'L')
    SeqLoc = Map.lSeq;
    RefSeqLoc = Map.lRefSeq;
end

%For each child seq, you want the RefSeq to be the Seq of the parent.
AncMap = renumberAncMap(AncMap); %Ensures AncMap numbering is relative to Tdata position. 
for j = 2:size(Tdata, 1)
    if AncMap(j, 2) == 0; continue; end %Leave root alone
    for c = 1:length(RefSeqLoc)
        Tdata(AncMap(j, 1), RefSeqLoc(c)) = Tdata(AncMap(j, 2), SeqLoc(c));
    end
end