%countSHM will count how many SHMs are present between the current sequence
%and the germline sequence of each group.
%
%  VDJdata = countSHM(VDJdata,NewHeader)

function VDJdata = countSHM(VDJdata,NewHeader)
getHeaderVar;

GrpNum = cell2mat(VDJdata(:,GrpNumLoc));
UnqGrpNum = unique(GrpNum);
for y = 1:length(UnqGrpNum)
    IdxLoc = find(GrpNum == UnqGrpNum(y));
    RefSeq = VDJdata{IdxLoc(1),RefSeqLoc};
    RefSeqXLoc = RefSeq == 'X';
    VMDNJ = cell2mat(VDJdata(IdxLoc(1),LengthLoc));
    
    if sum(VMDNJ) ~= length(RefSeq) || isempty(RefSeq)
        continue
    end
    
    for j = 1:length(IdxLoc)
        Seq = VDJdata{IdxLoc(j),SeqLoc};
        SeqXLoc = Seq == 'X';
        if length(Seq) ~= length(RefSeq); 
            warning('countSHM: Seq and RefSeq lengths do not match for entry # %d',IdxLoc(j));
            continue; 
        end
        MissLoc = Seq ~= RefSeq;
        MissLoc(SeqXLoc | RefSeqXLoc) = 0;
        
        Vmiss = sum(MissLoc(1:VMDNJ(1)));
        Mmiss = sum(MissLoc(VMDNJ(1)+1:sum(VMDNJ(1:2))));
        Dmiss = sum(MissLoc(sum(VMDNJ(1:2))+1:sum(VMDNJ(1:3))));
        Nmiss = sum(MissLoc(sum(VMDNJ(1:3))+1:sum(VMDNJ(1:4))));
        Jmiss = sum(MissLoc(sum(VMDNJ(1:4))+1:sum(VMDNJ(1:5))));
        
        if isempty(Vmiss); Vmiss = 0; end
        if isempty(Mmiss); Mmiss = 0; end
        if isempty(Dmiss); Dmiss = 0; end
        if isempty(Nmiss); Nmiss = 0; end
        if isempty(Jmiss); Jmiss = 0; end        
        
        VDJdata{IdxLoc(j),VmutLoc} = Vmiss;
        VDJdata{IdxLoc(j),MmutLoc} = Mmiss;
        VDJdata{IdxLoc(j),DmutLoc} = Dmiss;
        VDJdata{IdxLoc(j),NmutLoc} = Nmiss;
        VDJdata{IdxLoc(j),JmutLoc} = Jmiss;
    end
end