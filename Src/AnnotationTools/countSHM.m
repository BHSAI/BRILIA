%countSHM will count how many SHMs are present between the current sequence
%and the germline sequence of each group.
%
%  VDJdata = countSHM(VDJdata, Map)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    Map: structure map of BRILIA main data
%
%  OUTPUT
%    VDJdata: modified VDJdata where SHM count info are filled

function VDJdata = countSHM(VDJdata, Map)
GrpNum = cell2mat(VDJdata(:, Map.GrpNum));
UnqGrpNum = unique(GrpNum);
Idx = cell(length(UnqGrpNum), 1);
for y = 1:length(UnqGrpNum)
    Idx{y} = find(GrpNum == UnqGrpNum(y));
end

for c = 1:length(Map.Chain)
    Chain = lower(Map.Chain(c));
    SeqIdx  = Map.([Chain 'Seq']);
    RefSeqIdx  = Map.([Chain 'RefSeq']);
    LengthIdx  = Map.([Chain 'Length']);
    InvalidLoc = any(cellfun('isempty', VDJdata(:, [SeqIdx; RefSeqIdx; LengthIdx])), 2);
    
    VmutIdx = Map.([Chain 'Vmut']);
    NmutIdx = Map.([Chain 'Nmut']);
    JmutIdx = Map.([Chain 'Jmut']);
    if Chain == 'h'
        MmutIdx = Map.([Chain 'Mmut']);
        DmutIdx = Map.([Chain 'Dmut']);
    end
    
    for y = 1:length(UnqGrpNum)
        if InvalidLoc(Idx{y}(1)); continue; end

        RefSeq = VDJdata{Idx{y}(1), RefSeqIdx};
        SegLen = [VDJdata{Idx{y}(1), LengthIdx}];
        if sum(SegLen) ~= length(RefSeq); continue; end
        
        for j = 1:length(Idx{y})
            Seq = VDJdata{Idx{y}(j), SeqIdx};
            if length(Seq) ~= length(RefSeq)
                warning('%s: Seq and RefSeq lengths do not match for entry #%d', mfilename, Idx{y}(j));
                continue 
            end
            MissLoc = ~cmprSeqMEX(Seq, RefSeq, 'n');
            SumLen = cumsum(SegLen);
            
            VDJdata{Idx{y}(j), VmutIdx} = sum(MissLoc(1:SumLen(1)));
            VDJdata{Idx{y}(j), NmutIdx} = sum(MissLoc(SumLen(end-2)+1:SumLen(end-1)));
            VDJdata{Idx{y}(j), JmutIdx} = sum(MissLoc(SumLen(end-1)+1:SumLen(end)));
            if Chain == 'h'
                VDJdata{Idx{y}(j), MmutIdx} = sum(MissLoc(SumLen(1)+1:SumLen(2)));
                VDJdata{Idx{y}(j), DmutIdx} = sum(MissLoc(SumLen(2)+1:SumLen(3)));
            end
        end
    end
end