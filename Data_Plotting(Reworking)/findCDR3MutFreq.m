function FreqMat = findCDR3MutFreq(varargin)
if isempty(varargin)
    [VDJdata, VDJheader, ~, ~] = openSeqData;
else
    VDJdata = varargin{1};
    VDJheader = varargin{2};
end

H = getHeaderVar(VDJheader);

[Vmap, Dmap, Jmap] = getCurrentDatabase;

FreqMat = zeros(20,20);

%Separate 1-size clusters from >1-size clusters
H.GrpNumLoc = findHeader(VDJheader,'GroupNum');
GrpNum = cell2mat(VDJdata(:,H.GrpNumLoc));
GrpNumUnq = unique(GrpNum);
IdxMap = 1:size(GrpNum,1);

for y = 1:length(GrpNumUnq)
    if mod(y,10) == 0
        [y/length(GrpNumUnq)]
    end
    %Determine if group correction is possible
    GrpLoc = (GrpNumUnq(y) == GrpNum);
    IdxLoc = IdxMap(GrpLoc);
    Tdata = VDJdata(IdxLoc,:);    
    
    RefSeq = Tdata{1,H.RefSeqLoc};
    
    if sum(isnan(RefSeq))>0 || isempty(RefSeq); continue; end

    
    VMDNJ = cell2mat(VDJdata(IdxLoc(1),H.LengthLoc));
    Vdel = VDJdata{IdxLoc(1),H.DelLoc(1)};
    Jdel = VDJdata{IdxLoc(1),H.DelLoc(end)};    
    Vnum = VDJdata{IdxLoc(1),H.FamNumLoc(1)}(1);    
    Jnum = VDJdata{IdxLoc(1),H.FamNumLoc(end)}(1);    
    Cloc = VMDNJ(1) + Vdel - Vmap{Vnum,end} + 1;
    Wloc = sum(VMDNJ(1:4)) - Jdel + Jmap{Jnum,end} + 2; %Add 2 to capture the codon end.
    
    if Wloc > length(RefSeq)
        continue
    end
    
    CDR3Ref = nt2aa(RefSeq(Cloc:Wloc),'ACGTonly','false');
%     
%     [CDR3, ~] = findCDR3(Tdata(1,:),VDJheader,'RefSeq');
%     if isempty(CDR3{1})
%         continue
%     end
%     
% %     CDR3Ref = char(CDR3);   
%     CDR3Sam = cell(size(Tdata,1),1);
%     for w = 1:size(Tdata,1)
%         CDR3Sam{w} = nt2aa(Tdata{w,H.SeqLoc}(Cloc:Wloc));
%     end
%     CDR3Sam = char(CDR3Sam);
    
    CDR3Sam = Tdata(:,H.CDR3Loc(1));    
    DelThis = zeros(size(CDR3Sam,1),1)==1;
    for k = 1:size(CDR3Sam,1)    
        if ~isnan(CDR3Sam{k})
            if isempty(CDR3Sam{k}); 
                DelThis(k) = 1;
            end
        else
            DelThis(k) = 1;
        end
        
        if ~strcmpi(CDR3Sam{k}(1),'C') || ~strcmpi(CDR3Sam{k}(end),'W')
            DelThis(k) = 1;
        end
    end
    CDR3Sam(DelThis) = [];
    CDR3Sam = char(CDR3Sam);

    if isempty(CDR3Sam) || isempty(CDR3Ref)
        disp('No CDR3 Detected for group');
        continue
    end
    
    %Check for a * in the RefSeq.
    StarLoc = regexp(CDR3Ref,'\*');
    if ~isempty(StarLoc)
        CDR3SamCons = seqconsensus(CDR3Sam);
        for k = 1:length(StarLoc)
            CDR3Ref(StarLoc(k)) = CDR3SamCons(StarLoc(k));
        end
    end
    
    if size(CDR3Ref,2) ~= size(CDR3Sam,2)
        continue
    end

    [Pdis, Ploc] = calcPdist(CDR3Ref,CDR3Sam);
    
    %Find the maximum deviation sequence
    MaxDevLoc = Pdis == max(Pdis);
    CDR3Sam = CDR3Sam(MaxDevLoc,:);
    MutLoc = Ploc(MaxDevLoc,:);
    
    if sum(MutLoc) == 0; continue; end %No mutation found
    
    for q = 1:size(MutLoc,1)

        EvalRefMut = aa2int(CDR3Ref(MutLoc(q,:)));
        EvalRefMut(EvalRefMut > 20) = [];

        EvalSamMut = aa2int(CDR3Sam(q,MutLoc(q,:)));
        EvalSamMut(EvalSamMut > 20) = [];
    
        for k = 1:length(EvalSamMut)
            FreqMat(EvalSamMut(k),EvalRefMut(k)) = FreqMat(EvalSamMut(k),EvalRefMut(k)) + 1;
        end
    end
end
