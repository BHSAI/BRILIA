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
error('%s: this is obsolete. See clusterGeneByLineage', mfilename);

if nargin < 3
    GrpNumStart = 1;
end
Tdata(:, Map.ParNum) = {0};

RefSeqIdx = [Map.hRefSeq Map.lRefSeq];
RefSeqIdx(RefSeqIdx == 0) = [];
SeqIdx = [Map.hRefSeq Map.lRefSeq];
SeqIdx(SeqIdx == 0) = [];

%Determine the distances between sequences, adding H and L together
Chain = lower(Map.Chain);
PairDist = zeros(size(Tdata, 1));
for c = 1:length(Chain)
    Seq = Tdata(:, SeqIdx(c));
    [~, Motif, Mut, Penalty, ShmDist] = calcPairDistMEX(Seq);
    PairDist = PairDist + ShmDist;
    InvalidLoc = (Motif < 0 & Mut < 0) | Penalty >= 6; %Penalty of 6 means there's a double 2-consec miss. Penalty 9 means there's a triplet mismatch. Both are rare and should be discarded.
    PairDist(InvalidLoc) = Inf; %Prevents linking these
end
AncMap = calcAncMap(PairDist);

%Determine if you should break this up by G2C < P2C distance
Germ2SeqDist = zeros(size(AncMap, 1), 1);
for j = 1:size(AncMap, 1)
    for c = 1:length(Chain)
        [Ham, Motif, Mut, Penalty, ShmDist] = calcPairDistMEX(Tdata(AncMap(j, 1), [RefSeqIdx(c) SeqIdx(c)]));
        if (Motif(1, 2) < 0 && Mut(1, 2) < 0) || Penalty(1, 2) >= Ham(1, 2)
            Germ2SeqDist(j) = Inf;
        else
            Germ2SeqDist(j) = Germ2SeqDist(j) + ShmDist(1, 2);
        end
    end
end

AncMap = [AncMap findTreeClust(AncMap)];
for j = 1:max(AncMap(:, end))
    ClustLoc = AncMap(:, end) == j;
    if sum(ClustLoc) > 1
        RootSeqNum = AncMap(ClustLoc & Germ2SeqDist == min(Germ2SeqDist(ClustLoc)), 1); %Determine which one is the parent
        if length(RootSeqNum) > 1 %Break tie by higher template counts.
            TC = cell2mat(Tdata(RootSeqNum, Map.Template));
            RootSeqNum = RootSeqNum(TC == max(TC));
        end
        if length(RootSeqNum) > 1 %Break tie by one that's closest to all others
            PD = PairDist(RootSeqNum, :);
            PD(isinf(PD)) = 0;
            SumPD = sum(PD, 2);
            RootSeqNum = RootSeqNum(SumPD == min(SumPD));
        end
        AncMap(AncMap(:, 1) == RootSeqNum(1), 2) = 0; %If still tied, just get the first one
    end
end

%Now, make sure to break those that are closer to germline than neighbor
GermWinLoc = Germ2SeqDist < AncMap(:, 3);
AncMap(GermWinLoc, 2) = 0;
AncMap(GermWinLoc, 3) = Germ2SeqDist(GermWinLoc);
AncMap(:, 4) = findTreeClust(AncMap);

%Reroot each cluster according to the root
AncMap = sortrows(AncMap, [4 2 3 1]);
for j = 1:max(AncMap(:, end))
    ClustLoc = AncMap(:, end) == j;
    NumIdx = AncMap(ClustLoc, 1);
    if sum(ClustLoc) <= 2; continue; end
    RootedAncMap = calcRootedAncMap(PairDist(NumIdx, NumIdx));
    RootedAncMap(2:end, 2) = NumIdx(RootedAncMap(2:end, 2));
    RootedAncMap(:, 1) = NumIdx;
    AncMap(ClustLoc, 1:3) = RootedAncMap(:, 1:3);
end

%Rearrange TData, replacing RefSeq with Seq of parent, and giving it parent
%numbers. NOTE: in future releases, RefSeq will NOT be altered. 
Tdata = Tdata(AncMap(:, 1), :);
AncMap = renumberAncMap(AncMap);
for c = 1:length(Chain)
    Tdata(AncMap(:, 2) ~= 0, RefSeqIdx(c)) = Tdata(AncMap(AncMap(:, 2) ~= 0, 2), SeqIdx(c));
    Tdata(AncMap(:, 2) ~= 0, Map.ParNum) = Tdata(AncMap(AncMap(:, 2) ~= 0, 2), Map.SeqNum);
end
Tdata(:, Map.GrpNum) = num2cell(GrpNumStart - 1 + AncMap(:, end));

%Determine the number of children per each sequence
for k = 1:size(AncMap, 1)
    Tdata{k, Map.ChildCount} = length(findChild(AncMap, AncMap(k, 1)));
end
GrpNumStart = GrpNumStart + max(AncMap(:, end));

% The clustering of germline together was a bad idea in the sense that these tend to get grouped randomly, speicailly for shorter CDR3.
% It's better to separate them as suppose to try to join them. 
% 
% %The final step requires clustering germlines
% RootIdx = AncMap(AncMap(:, 2) == 0, 1);
% if numel(RefSeqIdx) == 1
%     RefSeq = Tdata(RootIdx, RefSeqIdx);
% else
%     RefSeq = cellfun(@(x, y) [x y], Tdata(RootIdx, RefSeqIdx(1)), Tdata(RootIdx, RefSeqIdx(2)));
% end
% 
% [RHam, RMotif, RMut, RPenalty, RShmDist] = calcPairDistMEX(RefSeq);
% RInvalidLoc = RMotif < 0 | RMut < 0 | RPenalty >= RHam; %This is stricter since you're clustering germlines!
% RShmDist(RInvalidLoc) = Inf;
% RAncMap = calcAncMap(RShmDist);
% RClustMap = findTreeClust(RAncMap);
% 
% %Renumber AncMap cluster number to be the same. Can have multiples 0's now
% GNum = 1;
% for k = 1:max(RClustMap)
%     LinkSeqNum = RootIdx(RClustMap == k);
%     [~, LinkSeqNumIdx] = intersect(AncMap(:, 1), LinkSeqNum, 'stable');
%     LinkGrpNumLoc = ismember(AncMap(:, 4), AncMap(LinkSeqNumIdx, 4));
%     AncMap(LinkGrpNumLoc, 4) = -GNum;
%     GNum = GNum + 1;
% end
% AncMap(:, 4) = abs(AncMap(:, 4));
