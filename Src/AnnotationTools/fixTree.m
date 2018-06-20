%fixTree will look through VDJdata tree structures,  and flip the root seq
%with first seq IF it doesn't make sense. One way to figure this out is to
%see if the 1st seq has MORE shm than another seq. Do this at the way end, 
%after D and N regions are fixed.
%
%  VDJdata = fixTree(VDJdata, VDJheader)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    Map: Structured map of BRILIA data types (getVDJmapper)
%
%  OUTPUT
%    VDJdata: modified VDJdata with cluster sequences reordered and
%      relinked via phylogeny lineage.

function VDJdata = fixTree(VDJdata, Map)

warning('%s: This is obsolete.', mfilename);
return

%Determine chain and needed locations
if strcmpi(Map.Chain, 'HL')
    SHMLoc = [Map.hVmut Map.hMmut Map.hDmut Map.hNmut Map.hJmut Map.lVmut Map.lNmut Map.lJmut];
    SeqLoc = [Map.hSeq Map.lSeq];
    RefSeqLoc = [Map.hRefSeq Map.lRefSeq];
elseif strcmpi(Map.Chain, 'H')
    SHMLoc = [Map.hVmut Map.hMmut Map.hDmut Map.hNmut Map.hJmut];
    SeqLoc = Map.hSeq;
    RefSeqLoc = Map.hRefSeq;
elseif strcmpi(Map.Chain, 'L')
    SHMLoc = [Map.lVmut Map.lNmut Map.lJmut];
    SeqLoc = Map.lSeq;
    RefSeqLoc = Map.lRefSeq;
end
ChildCountLoc = Map.ChildCount;
TemplateLoc = Map.Template;

%Begin fixing tree per cluster
GrpNum = cell2mat(VDJdata(:, Map.GrpNum));
UnqGrpNum = unique(GrpNum);
for y = 1:length(UnqGrpNum)
    %Extract necessary info
    IdxLoc = find(UnqGrpNum(y) == GrpNum);
    if length(IdxLoc) == 1; continue; end %Nothing to fix
    Tdata = VDJdata(IdxLoc, :);

    %Make sure you have the best root, lowest SHM total
    SHMCount = sum(cell2mat(Tdata(:, SHMLoc)), 2);
    MinIdx = find(SHMCount == min(SHMCount));
    if MinIdx(1) ~= 1 || length(MinIdx) > 1 %Need to double check
        TempCt = cell2mat(Tdata(MinIdx, TemplateLoc));        
        ChildCt = cell2mat(Tdata(MinIdx, ChildCountLoc));
        CompareMat = sortrows([MinIdx TempCt ChildCt], [-2 -3]); %Break tie by finding largest Template Ct or Child Ct
        MinIdx = CompareMat(1, 1);
    else
        MinIdx = MinIdx(1);
    end
    
    %Swap out the 1st seq of Tdata with the one that is closer to RefSeq
    if MinIdx ~= 1
        Tdata([1 MinIdx], :) = Tdata([MinIdx 1], :);
        for r = 1:length(RefSeqLoc) %Restore RefSeq, which should have been germline RefSeq
            Tdata([1 MinIdx], RefSeqLoc(r)) = Tdata([MinIdx 1], RefSeqLoc(r)); 
        end
    end

    %Proceed to rebuild tree. This generally improves lineages
    %Determine the distances between sequences
    PairDist = zeros(size(Tdata, 1));
    for c = 1:length(Map.Chain) %Add H and L chain distances
        [~, ShmPairDist] = calcPairDistMEX(Tdata(:, SeqLoc(c)));
        PairDist = PairDist + ShmPairDist;
    end
    
    %Calculate the new ancestral map
    AncMap = calcRootedAncMap(PairDist);
    Tdata = setAncRefSeq(Tdata, Map, AncMap); %This will remap RefSeq to Seq based on AncMap

    %Update the child count per parent
    ChildCt = zeros(size(AncMap, 1), 1);
    for k = 1:size(AncMap, 1)
        ChildCt(k) = length(findChild(AncMap, AncMap(k, 1)));
    end
    Tdata(:, ChildCountLoc) = num2cell(ChildCt);

    VDJdata(IdxLoc, :) = Tdata;
end
