%countSHM will count how many SHMs are present between the current sequence
%and the germline sequence of each group.
%
%  VDJdata = countSHM(VDJdata,VDJheader)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    VDJheader: main BRILIA header cell
%
%  OUTPUT
%    VDJdata: modified VDJdata where SHM count info are filled

function VDJdata = countSHM(VDJdata,Map)
if ~isstruct(Map) %Backward compatability
    Map = getVDJmapper(Map);
end
%Count SHMs per segment
GrpNum = cell2mat(VDJdata(:,Map.GrpNum));
UnqGrpNum = unique(GrpNum);
for y = 1:length(UnqGrpNum)
    IdxLoc = find(GrpNum == UnqGrpNum(y));
    
    for k = 1:length(Map.Chain)
        %Determine chain header locator
        if Map.Chain(k) == 'H'
            SeqLoc  = Map.hSeq;
            RefSeqLoc  = Map.hRefSeq;
            LengthLoc  = Map.hLength;
            VmutLoc = Map.hVmut;
            MmutLoc = Map.hMmut;
            DmutLoc = Map.hDmut;
            NmutLoc = Map.hNmut;
            JmutLoc = Map.hJmut;
        else
            SeqLoc  = Map.lSeq;
            RefSeqLoc  = Map.lRefSeq;
            LengthLoc  = Map.lLength;
            VmutLoc = Map.lVmut;
            NmutLoc = Map.lNmut;
            JmutLoc = Map.lJmut;
        end
        
        %Extract the RefSeq and segment info
        RefSeq = VDJdata{IdxLoc(1),RefSeqLoc};
        SegLen = cell2mat(VDJdata(IdxLoc(1),LengthLoc));
        
        %Make sure all information is provided for this data
        if isempty(RefSeq); continue; end
        if isempty(SegLen); continue; end        
        if sum(SegLen) ~= length(RefSeq); continue; end
        
        %Calculate SHM per each sequence in the group
        RefSeqXIdx = RefSeq == 'X';
        for j = 1:length(IdxLoc)
            %Extract remaining needed seq information
            Seq = VDJdata{IdxLoc(j),SeqLoc};
            
            %Make sure all information is provided for this data
            if isempty(Seq); continue; end
            if length(Seq) ~= length(RefSeq)
                warning('%s: Seq and RefSeq lengths do not match for entry #%d',mfilename,IdxLoc(j));
                continue; 
            end
            
            %Determine the sequence miss locations
            SeqXIdx = Seq == 'X';
            MissIdx = Seq ~= RefSeq;
            MissIdx(SeqXIdx | RefSeqXIdx) = 0;

            if Map.Chain(k) == 'H'
                try
                    Vmiss = sum(MissIdx(1:SegLen(1)));
                    Mmiss = sum(MissIdx(SegLen(1)+1:sum(SegLen(1:2))));
                    Dmiss = sum(MissIdx(sum(SegLen(1:2))+1:sum(SegLen(1:3))));
                    Nmiss = sum(MissIdx(sum(SegLen(1:3))+1:sum(SegLen(1:4))));
                    Jmiss = sum(MissIdx(sum(SegLen(1:4))+1:sum(SegLen(1:5))));
                catch
                    continue;
%                    save('debug_countSHM.mat')
                end

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
            else
                Vmiss = sum(MissIdx(1:SegLen(1)));
                Nmiss = sum(MissIdx(SegLen(1)+1:sum(SegLen(1:2))));
                Jmiss = sum(MissIdx(sum(SegLen(1:2))+1:sum(SegLen(1:3))));

                if isempty(Vmiss); Vmiss = 0; end
                if isempty(Nmiss); Nmiss = 0; end
                if isempty(Jmiss); Jmiss = 0; end        

                VDJdata{IdxLoc(j),VmutLoc} = Vmiss;
                VDJdata{IdxLoc(j),NmutLoc} = Nmiss;
                VDJdata{IdxLoc(j),JmutLoc} = Jmiss;
            end
        end
    end
end
