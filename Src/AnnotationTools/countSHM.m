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

function VDJdata = countSHM(VDJdata,VDJheader)
[H, L, Chain] = getAllHeaderVar(VDJheader);

%Count SHMs per segment
GrpNum = cell2mat(VDJdata(:,H.GrpNumLoc));
UnqGrpNum = unique(GrpNum);
for y = 1:length(UnqGrpNum)
    IdxLoc = find(GrpNum == UnqGrpNum(y));
    
    for k = 1:length(Chain)
        %Determine chain header locator
        if Chain(k) == 'H'
            B = H;
        else
            B = L;
        end
        
        %Extract the RefSeq and segment info
        RefSeq = VDJdata{IdxLoc(1),B.RefSeqLoc};
        SegLen = cell2mat(VDJdata(IdxLoc(1),B.LengthLoc));
        
        %Make sure all information is provided for this data
        if isempty(RefSeq); continue; end
        if isempty(SegLen); continue; end        
        if sum(SegLen) ~= length(RefSeq); continue; end
        
        %Calculate SHM per each sequence in the group
        RefSeqXIdx = RefSeq == 'X';
        for j = 1:length(IdxLoc)
            %Extract remaining needed seq information
            Seq = VDJdata{IdxLoc(j),B.SeqLoc};
            
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

            if Chain(k) == 'H'
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

                VDJdata{IdxLoc(j),B.VmutLoc} = Vmiss;
                VDJdata{IdxLoc(j),B.MmutLoc} = Mmiss;
                VDJdata{IdxLoc(j),B.DmutLoc} = Dmiss;
                VDJdata{IdxLoc(j),B.NmutLoc} = Nmiss;
                VDJdata{IdxLoc(j),B.JmutLoc} = Jmiss;
            else
                Vmiss = sum(MissIdx(1:SegLen(1)));
                Nmiss = sum(MissIdx(SegLen(1)+1:sum(SegLen(1:2))));
                Jmiss = sum(MissIdx(sum(SegLen(1:2))+1:sum(SegLen(1:3))));

                if isempty(Vmiss); Vmiss = 0; end
                if isempty(Nmiss); Nmiss = 0; end
                if isempty(Jmiss); Jmiss = 0; end        

                VDJdata{IdxLoc(j),B.VmutLoc} = Vmiss;
                VDJdata{IdxLoc(j),B.NmutLoc} = Nmiss;
                VDJdata{IdxLoc(j),B.JmutLoc} = Jmiss;
            end
        end
    end
end
