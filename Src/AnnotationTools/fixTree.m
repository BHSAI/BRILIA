%fixTree will look through VDJdata tree structures,  and flip the root seq
%with first seq IF it doesn't make sense. One way to figure this out is to
%see if the 1st seq has MORE shm than another seq. Do this at the way end, 
%after D and N regions are fixed.
%
%  VDJdata = fixTree(VDJdata, VDJheader)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    VDJheader: main BRILIA header cell
%
%  OUTPUT
%    VDJdata: modified VDJdata with cluster sequences reordered and
%      relinked via phylogeny lineage.

function VDJdata = fixTree(VDJdata, VDJheader)
%Determine chain and needed locations
[H, L, Chain] = getAllHeaderVar(VDJheader);
if strcmpi(Chain, 'HL')
    SHMLoc = [H.VmutLoc H.MmutLoc H.DmutLoc H.NmutLoc H.JmutLoc L.VmutLoc L.NmutLoc L.JmutLoc];
    SeqLoc = [H.SeqLoc L.SeqLoc];
    RefSeqLoc = [H.RefSeqLoc L.RefSeqLoc];
elseif strcmpi(Chain, 'H')
    SHMLoc = [H.VmutLoc H.MmutLoc H.DmutLoc H.NmutLoc H.JmutLoc];
    SeqLoc = H.SeqLoc;
    RefSeqLoc = H.RefSeqLoc;
elseif strcmpi(Chain, 'L')
    SHMLoc = [L.VmutLoc L.NmutLoc L.JmutLoc];
    SeqLoc = L.SeqLoc;
    RefSeqLoc = L.RefSeqLoc;
end
ChildCountLoc = max([H.ChildCountLoc L.ChildCountLoc]);
TemplateLoc = max([H.TemplateLoc L.TemplateLoc]);

%Begin fixing tree per cluster
GrpNum = cell2mat(VDJdata(:, H.GrpNumLoc));
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
    
    %Determine the maximum Seq Dist, which is used to set the class of PairDist
    MaxDist = 0;
    for c = 1:length(Chain)
        MaxDist = MaxDist + length(Tdata{1, SeqLoc(c)})^2 + length(Tdata{1, SeqLoc(c)})^2;
    end
    if MaxDist <= 2^8-1
        Class = 'uint8';
    elseif MaxDist <= 2^16-1
        Class = 'uint16';
    elseif MaxDist <= 2^32-1
        Class = 'uint32';
    else
        Class = 'double';
    end

    %Determine the distances between sequences
    PairDist = zeros(size(Tdata, 1), Class);
    for c = 1:length(Chain) %Add H and L chain distances
        PairDist = PairDist + calcPairDist(Tdata(:, SeqLoc(c)), 'shmham', Class);
    end
    
    %Calculate the new ancestral map
    AncMap = calcRootedAncMap(PairDist);
    Tdata = mapData2AncMap(Tdata, VDJheader, AncMap); %This will remap RefSeq to Seq based on AncMap

    %Update the child count per parent
    ChildCt = zeros(size(AncMap, 1), 1);
    for k = 1:size(AncMap, 1)
        ChildCt(k) = length(findChild(AncMap, AncMap(k, 1)));
    end
    Tdata(:, ChildCountLoc) = num2cell(ChildCt);

    VDJdata(IdxLoc, :) = Tdata;
end
