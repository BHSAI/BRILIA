%buildTreeLink will take a crude cluster of VDJdata (Tdata) and then use
%lineage trees to determine the fine clusters based on the SHM distance.
%This will return a rearranged Tdata cell matrix based on the new
%clustering, where the 1st member is always the sequence closest to the
%germline sequence.
%
%Procedure are as follows:
%1) Cluster based on distance first, cyclic dependencies allowed
%2) For each cluster, identify the root, defined as sequence with the
%   closest SHM distance to the annotated germline sequence.
%3) For each cluster, relink the tree from root and onward. 
%
%  Tdata = buildTreeLink(Tdata, Map)
%
%  INPUT
%    Tdata: subset of  the BRILIA data cell containing sequences that
%      belong in the same crude cluster (clusterGene.m) and have the
%      same sequence lengths (padtrimSeqGroup.m)
%    Map: map of the BRILIA header cells
%
%  OUTPUT
%    Tdata: rearranged Tdata such that the clusters are grouped together, 
%      with the 1st sequence of each cluster being closest to the germline
%      seq. For each seq, the field of "ParentNum" is filled with the
%      SeqNum of the immediate ancestor.

function [Tdata, GrpNumStart] = buildTreeLink(Tdata, Map, GrpNumStart)
if nargin < 3
    GrpNumStart = 1;
end
Tdata(:, Map.ParNum) = {0};

%Determine the distances between sequences, adding H and L together
Chain = lower(Map.Chain);
PairDist = zeros(size(Tdata, 1));
for c = 1:length(Chain)
    Seq = Tdata(:, Map.([Chain(c) 'Seq']));
    [Ham, Motif, Mut, Penalty] = calcSeqShmMEX(Seq);
    InvalidLoc = (Motif < 0 & Mut < 0) | Penalty./Ham > 0.5;
    PairDist = PairDist + (Ham -(Motif + Mut - Penalty)/4);
    PairDist(InvalidLoc) = Inf; %Prevents linking these
end

AncMap = calcAncMap(PairDist);
AncMap = [AncMap findTreeClust(AncMap)];

for j = 1:max(AncMap(:, end))
    ClustLoc = AncMap(:, end) == j;
    if sum(ClustLoc) <= 1
        continue
    end
    
    %Determine which one is the parent
    RootSeqNum = AncMap(ClustLoc, 1);
    Germ2SeqDist = zeros(length(RootSeqNum), 1);
    for k = 1:length(RootSeqNum)
        for c = 1:length(Chain)
            [Ham, Motif, Mut, Penalty] = calcSeqShmMEX(Tdata(RootSeqNum(k), [Map.([Chain(c) 'Seq']) Map.([Chain(c) 'RefSeq'])]));
            Germ2SeqDist(k) = Germ2SeqDist(k) + Ham(1, 2) - (Motif(1, 2) + Mut(1, 2) - Penalty(1, 2))/4;
        end
    end
    RootSeqNum = RootSeqNum(Germ2SeqDist == min(Germ2SeqDist));
        
    %Break tie by higher template counts
    if length(RootSeqNum) > 1
        TC = cell2mat(Tdata(RootSeqNum, Map.Template));
        RootSeqNum = RootSeqNum(TC == max(TC));
    end

    %Break tie by one that's closest to all others
    if length(RootSeqNum) > 1
        PD = PairDist(RootSeqNum, :);
        PD(isinf(PD)) = 0;
        SumPD = sum(PD, 2);
        RootSeqNum = RootSeqNum(SumPD == min(SumPD));
    end

    AncMap(AncMap(:, 1) == RootSeqNum(1), 2) = 0; %If still tied, just get the first one
end

AncMap = sortrows(AncMap, [4 2]);

%Reroot each cluster according to the root
for j = 1:max(AncMap(:, end))
    ClustLoc = AncMap(:, end) == j;
    NumIdx = AncMap(ClustLoc, 1);
    if sum(ClustLoc) <= 2; continue; end
    RootedAncMap = calcRootedAncMap(PairDist(NumIdx, NumIdx));
    RootedAncMap(2:end, 2) = NumIdx(RootedAncMap(2:end, 2));
    RootedAncMap(:, 1) = NumIdx;
    AncMap(ClustLoc, 1:3) = RootedAncMap(:, 1:3);
end

%For each root, see if it is closer to its germline or another cluster
%member.



%Rearrange TData, replacing RefSeq with Seq of parent, and giving it parent
%numbers. NOTE: in future releases, RefSeq will NOT be switched off. 
Tdata = Tdata(AncMap(:, 1), :);
AncMap = renumberAncMap(AncMap);
for c = 1:length(Chain)
    Tdata(AncMap(:, 2) ~= 0, Map.([Chain(c) 'RefSeq'])) = Tdata(AncMap(AncMap(:, 2) ~= 0, 2), Map.([Chain(c) 'Seq']));
    Tdata(AncMap(:, 2) ~= 0, Map.ParNum) = Tdata(AncMap(AncMap(:, 2) ~= 0, 2), Map.SeqNum);
end
Tdata(:, Map.GrpNum) = num2cell(GrpNumStart - 1 + AncMap(:, end));

%Determine the number of children per each sequence
for k = 1:size(AncMap, 1)
    Tdata{k, Map.ChildCount} = length(findChild(AncMap, AncMap(k, 1)));
end
GrpNumStart = GrpNumStart + max(AncMap(:, end));