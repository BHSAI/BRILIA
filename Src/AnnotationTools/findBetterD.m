%findBetterD will attempt to find a better D gene match per group by trying
%to increase the VDJ alignment scores and Nscores.
%
%Nscore = (2*Ptdt-1)L^2 where Ptdt is the TDT score and L is the length of
%the N region.
%
%NOTE: This Nscore calculation differs from BRILIA paper publication
%equation. It was changed to this one to reflect the difference in
%probability of random NT vs TDT NT addition. The original equation, for
%instance, will yield high alignment score for long sequences with random
%nts. This one will yield negative values if that happens.
%
%  VDJdata = findBetterD(VDJdata, VDJheader, DB)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    VDJheader: main BRILIA header cell
%    DB: Gene database structure (getCurrentDatabase.m)
%
%  OUTPUT
%    VDJdata: modified VDJdata where D genes are recomputed that maximizes
%      the VDJ alignment scores and N region scores.
%
%  See also calcTDTscore

function VDJdata = findBetterD(VDJdata, Map, DB)
if strcmpi(Map.Chain, 'L')
    return
end

%Look for better D
GrpNum = cell2mat(VDJdata(:, Map.GrpNum));
UnqGrpNum = unique(GrpNum);
UpdateIdx = zeros(size(VDJdata, 1), 1, 'logical');
for y = 1:length(UnqGrpNum)
    if ~mod(y, 200)
        showStatus(sprintf('  Refining D gene  %d / %d.', y, length(UnqGrpNum)));
    end
    
    %Identify group, if any. For single, check all mutations.
    IdxLoc = find(UnqGrpNum(y) == GrpNum);
    Tdata = VDJdata(IdxLoc, :);
    VMDNJ = cell2mat(Tdata(1, Map.hLength));

    %Extract necessary V informations
    Vname = Tdata{1, Map.hGeneName(1)};
    VdelCur = Tdata{1, Map.hDel(1)};
    VmapNum = Tdata{1, Map.hGeneNum(1)};
    VallowedDel = DB.Vmap{VmapNum(1), end} - 3; %Correct -3 as deletion length is AFTER C codon.
    if VallowedDel <= 0; VallowedDel = 25; end

    %Extract necessary J informations
    Jname = Tdata{1, Map.hGeneName(3)};
    JdelCur = Tdata{1, Map.hDel(end)};
    JmapNum = Tdata{1, Map.hGeneNum(end)};
    JallowedDel = DB.Jmap{JmapNum(1), end} - 1; %Correct -1 as deletion length is BEFORE F/W codon.
    if JallowedDel < 0; JallowedDel = 25; end

    %Find the mismatched nts with respect to 1st seq of cluster only.
    RefSeq = Tdata{1, Map.hRefSeq};
    XlocRef = RefSeq == 'X';
    ConsMissCt = zeros(size(RefSeq));
    ErrorDetected = 0;
    for k = 1:size(Tdata, 1)
        Seq = Tdata{k, Map.hSeq};
        if length(Seq) ~= length(RefSeq)
            ErrorDetected = 1;
            break;
        end
        XlocSeq = Seq == 'X';
        MissLoc = ~(RefSeq == Seq | XlocSeq | XlocRef);
        ConsMissCt = ConsMissCt + MissLoc;
    end    
    if ErrorDetected
        warning('%s: Erroroenous annotation detected.', mfilename);
        continue;
    end

    %Find consensus misses that are >= MaxMiss as seen in the V segment,
    %which is the most accurate segment
    VconsSeg = VMDNJ(1) + VdelCur - VallowedDel;
    MaxMiss = max(ConsMissCt(1:VconsSeg));
    if size(Tdata, 1) == 1 || isempty(MaxMiss)
        ConsMiss = ConsMissCt >= 1; %Must be atleast 1
    else
        ConsMiss = ConsMissCt >= MaxMiss;
    end

    %Determine location of mismatches after 104C
    VcutLen = sort(VMDNJ(1) - find(ConsMiss(1:VMDNJ(1)) == 1) + 1);
    VcutLen(VcutLen > VMDNJ(1)) = []; %Just in case there is no cons mismatch
    VtotDel = VcutLen + VdelCur; %How much of the reference gene must be cut to achieve this cutLoc.
    VcutLen(VtotDel > VallowedDel) = []; %Remaining allowed cut len.
    if isempty(VcutLen)
        VcutLen = 0;
    end

    %Determine location of mismatches right before 118W. 
    JcutLen = find(ConsMiss(end-VMDNJ(end)+1:end) == 1);
    JtotDel = JcutLen + JdelCur;
    JcutLen(JtotDel > JallowedDel) = [];
    if isempty(JcutLen)
        JcutLen = 0;
    end

    if max(VcutLen) > 0 || max(JcutLen) > 0 %Change probably needed

        %Make sure cutLen is not empty, and has a 0.
        if VcutLen(1) ~= 0
            VcutLen = cat(2, 0, VcutLen);
        end
        if JcutLen(1) ~= 0 
            JcutLen = cat(2, 0, JcutLen);
        end

        %Calculate the various D alignment result for all combination of Vcut    
        CompareMat = zeros(length(JcutLen)*length(VcutLen), 7); %[Vcut Jcut Dscore Vscore Jscore NvdScore NdjScore]
        Dmatch = cell(length(JcutLen)*length(VcutLen), 6);
        MissRate = MaxMiss/VconsSeg;

        q = 1;
        for v = 1:length(VcutLen)
            for j = 1:length(JcutLen)
                TestDseg = RefSeq(VMDNJ(1)-VcutLen(v)+1:end-VMDNJ(end)+JcutLen(j));
                AllowedMiss = ceil(MissRate * length(TestDseg));

                Dmatch(q, :) = findGeneMatch(TestDseg, DB.Dmap, 'D', AllowedMiss);
                CompareMat(q, 1) = VcutLen(v);
                CompareMat(q, 2) = JcutLen(j);
                CompareMat(q, 3) = Dmatch{q, 5}(2); %Scoring results

                %Calculate the new V score
                VconsMatch = ~ConsMiss(1:VMDNJ(1)-VcutLen(v));
                CompareMat(q, 4) = calcAlignScoreMEX(VconsMatch);                

                %Calculate the new J score
                JconsMatch = ~ConsMiss(end-VMDNJ(end)+1+JcutLen(j):end);
                CompareMat(q, 5) = calcAlignScoreMEX(JconsMatch);

                %Get the new Nvd D Ndj lengths
                Mlen = Dmatch{q, 4}(1);
                Dlen = Dmatch{q, 4}(2);
                Nlen = Dmatch{q, 4}(3);

                %Calculate TDT score for Nvd
                Mseq = '';
                if Mlen >= 1
                    Mseq = RefSeq(VMDNJ(1)-VcutLen(v)+1:VMDNJ(1)-VcutLen(v)+Mlen);
                end
                MtdtScore = calcTDTscore(Mseq);
                if isempty(MtdtScore); MtdtScore = 0; end

                %Calculate TDT score for Ndj
                Nseq = '';
                if Nlen >= 1
                    Nseq = RefSeq(VMDNJ(1)-VcutLen(v)+Mlen+Dlen+1:VMDNJ(1)-VcutLen(v)+Mlen+Dlen+Nlen);
                end
                NtdtScore = calcTDTscore(Nseq);
                if isempty(NtdtScore); NtdtScore = 0; end

                %Calculate Nscores for Nvd and Ndj
                CompareMat(q, 6) = (2*MtdtScore - 1)*length(Mseq)^2;
                CompareMat(q, 7) = (2*NtdtScore - 1)*length(Nseq)^2;

                q = q+1;
            end
        end

        %Determine maximum alignment score for D
        TotScore = sum(CompareMat(:, 3:7), 2);
        BestD = TotScore == max(TotScore);  
        if sum(BestD) > 1 %Break ties by looking at the D end deletion counts
            D5D3del = cell2mat(Dmatch(:, 3));
            MaxDels = max(D5D3del(:, [1 3]), [], 2);
            BestD = BestD & (MaxDels == min(MaxDels(BestD)));
        end    
        BestMatch = find(BestD == 1);
        BestMatch = BestMatch(1); %If still tie, take 1st one only.
        if BestMatch == 1; continue; end %no changes needed after all

        %Update the necessary informations for the Tdata
        Dmatch = Dmatch(BestMatch, :);
        VnewDel = CompareMat(BestMatch, 1); %Nts to trim from V portion
        JnewDel = CompareMat(BestMatch, 2); %Nts to trim from J portion
        VMDNJnew = [VMDNJ(1)-VnewDel  Dmatch{1, 4}  VMDNJ(end)-JnewDel];

        %Attempt a Vmatch realignment,  using leftover as M's.
        if VMDNJnew(2) > 0 
            Vnt = RefSeq(1:sum(VMDNJnew(1:2)));
            AllowedMiss = ceil(MissRate * length(Vnt));

            CDR3start = Tdata{1, Map.hCDR3(3)};
            if isempty(CDR3start); CDR3start = 0; end
            Vmatch = findGeneMatch(Vnt, DB.Vmap, 'V', AllowedMiss, CDR3start); %Redo for all V's

            VMDNJnew(1) = sum(Vmatch{4}(1:2));
            VMDNJnew(2) = Vmatch{4}(3);
            VnewDel = Vmatch{3}(1, 3) - VdelCur; %Nts to trim/add from V. Subtract VdelCur, since you'll add it back.
            VmapNum = Vmatch{1};
            Vname = Vmatch{2};
        end

        %Attempt a Jmatch realignment, using leftover as N's.
        if VMDNJnew(4) > 0 
            Jnt = RefSeq(sum(VMDNJnew(1:3))+1:end);
            AllowedMiss = ceil(MissRate * length(Jnt));

            CDR3end = Tdata{1, Map.hCDR3(4)} - sum(VMDNJnew(1:3)); %need to adjust for Jnt being smaller than first Seq
            if isempty(CDR3end); CDR3end = 0; end
            Jmatch = findGeneMatch(Jnt, DB.Jmap, 'J', AllowedMiss, CDR3end); %Redo for all J's.

            VMDNJnew(4) = Jmatch{4}(1);
            VMDNJnew(5) = sum(Jmatch{4}(2:3));        
            JnewDel = Jmatch{3}(1, 1) - JdelCur; %Nts to trim/add from J.  Subtract JdelCur, since you'll add it back.
            JmapNum = Jmatch{1};
            Jname = Jmatch{2};
        end

        %Make sure the new D isn't returning invalid VMDNJ lengths
        if min(VMDNJnew([1 3 5])) <= 0 %Try not to lose the VDJ.
            warning('%s: Skipping D fix due to loss of V, D, or J.', mfilename);
            continue;
        end 
        if min(VMDNJnew([2 4])) < 0 %Don't fix errors
            warning('%s: Skipping D fix due to negative Nvd or Ndj.', mfilename);
            continue; 
        end 

        %Update VDJdata and mark which ones need RefSeq / SHM updates
        VDDJdels = [(VdelCur+VnewDel)  Dmatch{1, 3}(1, 1)  Dmatch{1, 3}(1, 3)  (JdelCur+JnewDel)];
        Tdata(:, Map.hDel) = repmat(num2cell(VDDJdels), size(Tdata, 1), 1);
        Tdata(:, Map.hLength) = repmat(num2cell(VMDNJnew), size(Tdata, 1), 1);
        Tdata(:, Map.hGeneNum) = repmat({VmapNum Dmatch{1, 1} JmapNum}, size(Tdata, 1), 1);
        Tdata(:, Map.hGeneName) = repmat({Vname Dmatch{1, 2} Jname}, size(Tdata, 1), 1);
        VDJdata(IdxLoc, :) = Tdata;
        UpdateIdx(IdxLoc) = 1;
    end
end

%Update those that have changed
VDJdata(UpdateIdx, :) = buildRefSeq(VDJdata(UpdateIdx, :), Map, DB, 'H', 'germline', 'first'); %must do first seq of all cluster
VDJdata(UpdateIdx, :) = updateVDJdata(VDJdata(UpdateIdx, :), Map, DB);