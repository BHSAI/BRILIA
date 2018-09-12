%clusterByLineage will cluster a group according to lineage linkability.
%This uses hamming distance, motif, substitution direction, and consecutive
%mismatch penalty to determine if two clones belong to the same clonotype. 
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    Map: header map of BRILIA data cells
%
%  OUTPUT
%    VDJdata: Grp numbers are reassigned based on the lineage linkability
%
%  NOTE
%    This algorithm is designed to work WITHOUT any hamming distance cutoff
%    criteria. However, we do impose a 10% hamming distance cutoff for the
%    "PENALTY" for consectuve mismatch within the CDR3 length. That is, if
%    there are consecutive mismatch, it'll generate a Mi^2-1 penalty. If
%    this penalty is > 10% of the CDR3 nt length, then it'll prevent
%    linking under the same clonotype.

function VDJdata = clusterByLineage(VDJdata, Map, varargin)
if isempty(VDJdata); return; end

UseSHMHAM = startsWith('shmham', varargin, 'ignorecase', true);
KeepCyclic = startsWith('cyclic', varargin, 'ignorecase', true);
IsSpliced = iscell(VDJdata{1});
if ~IsSpliced
    VDJdata = spliceData(VDJdata, Map); %Has to be done to ensure proper group numbers
end
parfor y = 1:length(VDJdata)
    VDJdata{y} = clusterByLineagePerGroup(VDJdata{y}, Map, UseSHMHAM, KeepCyclic);
end
VDJdata = joinData(VDJdata, Map); %Have to join to ensure cluster-based splicing w/ correct GrpNum
if IsSpliced
    VDJdata = spliceData(VDJdata, Map); %Re-splice to ensure cluster-based splicing
end

function Tdata = clusterByLineagePerGroup(Tdata, Map, UseSHMHAM, KeepCyclic)
if isempty(Tdata) || size(Tdata, 1) == 1; return; end
LongDist = 1E5; %Distance to use to prevent linking. Can't be Inf to other reasons.
Tdata = padtrimSeqGroup(Tdata, Map, 'cdr3length', 'max', 'Seq');
if isempty(Tdata); return; end

InvalidLoc = zeros(size(Tdata, 1), 'logical');
PairDist = zeros(size(Tdata, 1));
for c = 1:length(Map.Chain)
    C = lower(Map.Chain(c));
    Seq  = Tdata(:, Map.([C 'Seq']));
    [CDR3Len, CDR3S, CDR3E] = deal(Tdata{1, Map.([C 'CDR3'])(2:4)});
    CDR3 = cellfun(@(x) x(CDR3S:CDR3E), Seq, 'un', 0);
    [HamDist1, Motif1, Mut1, Penalty1, ShmDist1] = calcPairDistMEX(Seq);
    [~, Motif2, Mut2, Penalty2] = calcPairDistMEX(CDR3);
    InvalidLoc1 = ((Motif1 < 0) & (Mut1 < 0)) | (Penalty1/length(Seq) >= 0.1);
    InvalidLoc2 = ((Motif2 < 0) & (Mut2 < 0)) | (Penalty2/CDR3Len     >= 0.1);
    InvalidLoc = InvalidLoc | InvalidLoc1 | InvalidLoc2;
    if UseSHMHAM
        PairDist = PairDist + ShmDist1;
    else
        PairDist = PairDist + HamDist1;
    end
end

%If there is a 0 distance (regardless of Ham or ShmHam dist), this means the sequences are same. Group them up.
[X, Y] = find(PairDist == 0);
DiagLoc = X == Y;
X = unique(X(~DiagLoc));
KeepLoc = repelem(true, 1, size(PairDist, 2));
if ~isempty(X) %Need to modify
    KeepLoc(X(2:end)) = 0;
    Tdata{X(1), Map.Template} = sum(cell2mat(Tdata(X, Map.Template)));
    Tdata = Tdata(KeepLoc, :);
    PairDist = PairDist(KeepLoc, KeepLoc);
    InvalidLoc = InvalidLoc(KeepLoc, KeepLoc);
end

if KeepCyclic
    AncMap = calcAncMap(PairDist);
    AncMap(:, 4) = findTreeClust(AncMap);
else
    PairDist(InvalidLoc) = LongDist; %make it high, but not inf;

    %Determine if you should break this up by G2C < P2C distance
    G2CDist = zeros(size(Tdata, 1), 1);
    for j = 1:size(Tdata, 1)
        for c = 1:length(Map.Chain)
            C = lower(Map.Chain(c));
            [HamDist, ~, ~, ~, ShmDist] = calcPairDistMEX(Tdata(j, [Map.([C 'RefSeq']) Map.([C 'Seq'])]));
            if UseSHMHAM
                G2CDist(j) = G2CDist(j) + ShmDist(1, 2);
            else
                G2CDist(j) = G2CDist(j) + HamDist(1, 2);
            end
        end
    end

    AncMap = calcAncMap(PairDist);

    %Makesure to break up anything with longer distance into separate clonotypes
    CutLoc = AncMap(:, 3) >= LongDist;
    AncMap(CutLoc, 2) = 0;
    AncMap(CutLoc, 3) = G2CDist(CutLoc);
    AncMap(:, 4) = findTreeClust(AncMap);

    for j = 1:max(AncMap(:, end))
        ClustLoc = AncMap(:, end) == j;
        if any(ClustLoc)
            RootSeqNum = AncMap(ClustLoc & G2CDist == min(G2CDist(ClustLoc)), 1); %Determine which one is the parent
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
    GermWinLoc = G2CDist < AncMap(:, 3);
    AncMap(GermWinLoc, 2) = 0;
    AncMap(GermWinLoc, 3) = G2CDist(GermWinLoc);
    AncMap(:, 4) = findTreeClust(AncMap);

    %Reroot each cluster according to the root
    AncMap = sortrows(AncMap, [4 2 3 1]);
    for j = 1:max(AncMap(:, end))
        ClustLoc = AncMap(:, end) == j;
        NumIdx = AncMap(ClustLoc, 1);
        RootedAncMap = calcRootedAncMap(PairDist(NumIdx, NumIdx));
        RootedAncMap(2:end, 2) = NumIdx(RootedAncMap(2:end, 2));
        RootedAncMap(:, 1) = NumIdx;
        AncMap(ClustLoc, 1:3) = RootedAncMap(:, 1:3);
    end

    if any(findTreeCycle(AncMap))
        error('%s: something is wrong');
        warning('%s: Found a cyclic dependency in group number %d.', mfilename, UnqGrpNum);
    end
    
    Tdata = Tdata(AncMap(:, 1), :);
    AncMap = renumberAncMap(AncMap);

end
    
%Reassign the RefSeq as the Seq of the parent. NOTE: in future releases, RefSeq will NOT be altered. 
for c = 1:length(Map.Chain)
    C = lower(Map.Chain(c));
    Tdata(AncMap(:, 2) ~= 0, Map.([C 'RefSeq'])) = Tdata(AncMap(AncMap(:, 2) ~= 0, 2), Map.([C 'Seq']));
end
Tdata(AncMap(:, 2) ~= 0, Map.ParNum) = Tdata(AncMap(AncMap(:, 2) ~= 0, 2), Map.SeqNum);

%Determine the number of children per each sequence
for k = 1:size(AncMap, 1)
    Tdata{k, Map.ChildCount} = length(findChild(AncMap, AncMap(k, 1)));
end

Tdata(:, Map.GrpNum) = num2cell(AncMap(:, 4));
