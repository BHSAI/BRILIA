%findBetterD will attempt to find a better D gene match per group by
%exploiting the locations of consensus mismatches near the V and J edges.
%
% VDJdata = findBetterD(VDJdata,NewHeader,Vmap,Dmap,Jmap)
%
% 1) Identify "consensus" mismatch using majority rule, inclusive of 50%
% split.
%
% 2) Always trim at >= 3 consec mismatch and happens between 104C and 118W
% region.
%     EX: |||xxx||  --> oooooo|| if we are trimming left side, or
%                   --> |||ooooo if we are trimming right side.
% 
% 3) Always cut on concensus mismatches with the same NT for group >= 3
%     EX: 'ACGTGgGT'   will become   'ACGTG ggt'
%         'ACcTGgGT'   will become   'ACcTG ggt'
%         'AcGTGgGT'   will become   'AcGTG ggt'
%         'ACaTGgGT'   will become   'ACaTG ggt'
%         'ACGaGgGT'   will become   'ACGaG ggt'
%
% 4) Decide to trim for <= 2 consec mismatch based on how well it
% improves/worsen alignment score, for majority mismatch rule.
%     EX: |x||||||  --> oo|||||| 
%     EX: ||x|||||  --> ||x||||| no trim
%     EX: ||xx||||  --> oooo||||
%     EX: |x|x||||  --> oooo||||
%     EX: |||xx|||  --> |||xx||| no trim


% For TroubleShooting
% EvalRange = 21;
% [VDJdata, NewHeader] = openSeqData;
% VDJdata = VDJdata(EvalRange,:);
% VDJdata = reformatAlignment(VDJdata,3);

function VDJdata = findBetterD(VDJdata,NewHeader,varargin)
%Extract the VDJ database
if length(varargin) == 3
    Vmap = varargin{1};
    Dmap = varargin{2};
    Jmap = varargin{3};
else
    [Vmap, Dmap, Jmap] = getCurrentDatabase;
end

getHeaderVar;

%Determine V and J pad lengths, based on allowed deletion counts. Helps
%with overhang matches by convolveSeq.
VpadCt = max(cell2mat(Vmap(:,end)));
JpadCt = max(cell2mat(Jmap(:,end)));
Vpad = repmat('*',1,VpadCt);
Jpad = repmat('*',1,JpadCt);

%Look for better D
GrpNum = cell2mat(VDJdata(:,GrpNumLoc));
GrpNumUnq = unique(GrpNum);
IdxMap = 1:size(GrpNum,1);
for y = 1:length(GrpNumUnq)
    if mod(y,20) == 0
        [y/length(GrpNumUnq)]
    end
    %Identify group, if any. For single, check all mutations.
    GrpLoc = (GrpNumUnq(y) == GrpNum);
    IdxLoc = IdxMap(GrpLoc);
    Tdata = VDJdata(GrpLoc,:);      %Temporarily pullout group data and necesary informations
    
    %Extract necessary informations
    [VMDNJ, ~, ~] = unique(cell2mat(Tdata(:,LengthLoc)),'rows');
    
    VdelCur = Tdata{1,DelLoc(1)};
    VmapNum = Tdata{1,FamNumLoc(1)};
    Vname = Tdata{1,FamLoc(1)};
    VallowedDel = Vmap{VmapNum(1),end} - 3; %Correct -3 as deletion length is AFTER C codon.
    if VallowedDel < 0; VallowedDel = 25; end    
    
    JdelCur = Tdata{1,DelLoc(end)};
    JmapNum = Tdata{1,FamNumLoc(end)};
    Jname = Tdata{1,FamLoc(3)};
    JallowedDel = Jmap{JmapNum(1),end} - 1; %Correct -1 as deletion length is BEFORE F/W codon.
    if JallowedDel < 0; JallowedDel = 25; end
    
    %Find the location of ANY mismatches within the group
    FirstSeq = Tdata{1,SeqLoc}; %clusterGene will have placed this as the 1st seqeuence of group as closest to Ref Seq.
    Classifier = char(Tdata(:,FormClassLoc(2)));
    ConsMiss = zeros(1,size(Classifier,2));
    for k = 1:size(Classifier,1)
        ConsMiss = ConsMiss + isstrprop(Classifier(k,:),'lower');
    end
    
    %Find max # of consensus mismatch possible based on the V before 104C
    if size(Classifier,1) > 1
        ConsMissVb4C = ConsMiss(1:VMDNJ(1)+VdelCur-VallowedDel);
        if isempty(ConsMissVb4C); ConsMissVb4C = 1; end
        ConsNum = max(ConsMissVb4C);
        if ConsNum == 0; ConsNum =1; end %Must be at least 1.
        ConsMiss = ConsMiss >= ConsNum; %Only consensus mismatches >= # found in V gene are considered
    end

    %Determine location of mismatches after 104C
    VcutLen = sort(VMDNJ(1) - find(ConsMiss(1:VMDNJ(1)) == 1) + 1);
    VcutLen(VcutLen > VMDNJ(1)) = []; %Just in case there is no cons mismatch
    VtotDel = VcutLen + VdelCur; %How much of the reference gene must be cut to achieve this cutLoc.
    VcutLen(VtotDel > VallowedDel) = []; %Remaining allowed cut len.
    
    %Determine location of mismatches right before 118W. 
    JcutLen = find(ConsMiss(end-VMDNJ(end)+1:end) == 1);
    JtotDel = JcutLen + JdelCur; %How much of the reference gene must be cut to achieve this cutLoc.
    %JcutLen(JtotDel > JallowedDel) = []; %Remaining allowed cut len.
    
    if ~isempty(VcutLen) || ~isempty(JcutLen) %Change probably needed

        %Make sure cutLen is not empty, and has a 0.
        if isempty(VcutLen); 
            VcutLen = 0;
        else
            VcutLen = [0 VcutLen];
        end
        if isempty(JcutLen); 
            JcutLen = 0; 
        else
            JcutLen = [0 JcutLen];
        end
        
        MissRate = Tdata{1,SHMLoc(1)}/VMDNJ(1);
            
        %Calculate the various D alignment result for all combination of Vcut    
        CompareMat = zeros(length(JcutLen)*length(VcutLen),7); %[Vcut Jcut Dscore Vscore Jscore N2length N1length]
        Dmatch = cell(length(JcutLen)*length(VcutLen),6);
        q = 1;
        for v = 1:length(VcutLen)
            for j = 1:length(JcutLen)
                NTeval = FirstSeq(:,VMDNJ(1)-VcutLen(v)+1:end-VMDNJ(end)+JcutLen(j));
                AllowedMiss = ceil(MissRate * length(NTeval));

                Dmatch(q,:) = findGeneMatch(NTeval,Dmap,'D',AllowedMiss);
                CompareMat(q,1) = VcutLen(v);
                CompareMat(q,2) = JcutLen(j);
                CompareMat(q,3) = Dmatch{q,5}(1,2); %Scoring results
                
                %Calculate the new V score
                VconsMatch = (ConsMiss(1:VMDNJ(1)-VcutLen(v)) == 0);
                CompareMat(q,4) = calcAlignScore(VconsMatch);                
                
                %Calculate the new J score
                JconsMatch = (ConsMiss(end-VMDNJ(end)+1+JcutLen(j):end) == 0);
                CompareMat(q,5) = calcAlignScore(JconsMatch);
                
                %Save the N lengths
                Dlen = Dmatch{q,4}(2);
                N2len = Dmatch{q,4}(1);
                if N2len > 1
                    N2seq = FirstSeq(VMDNJ(1)-VcutLen(v)+1:VMDNJ(1)-VcutLen(v)+N2len);
                else
                    N2seq = '';
                end
                
                N1len = Dmatch{q,4}(3);
                if N1len > 1
                    N1seq = FirstSeq(VMDNJ(1)-VcutLen(v)+N2len+Dlen+1:VMDNJ(1)-VcutLen(v)+N2len+Dlen+N1len);
                else
                    N1seq = '';
                end
                                
                %--
%                 Pn2 = calcTDTscore(N2seq);
%                 Ln2 = length(N2seq);
%                 Pn1 = calcTDTscore(N1seq);
%                 Ln1 = length(N1seq);
%                 
%                 CompareMat(q,6) = Ln2^2*(Pn2^2-(1-Pn2)^2); %))(1-Pn2)((2*Pn2-1)*length(N2seq))^2; %Dmatch{q,4}(1); %N2 region length
%                 CompareMat(q,7) = Ln1^2*(Pn1^2-(1-Pn1)^2); %((2*Pn1-1)*length(N1seq))^2; %Dmatch{q,4}(3); %N1 region length
%                 
%                 
                %--
%                 CompareMat(q,6) = (calcTDTscore(N2seq)*length(N2seq))^2; %Dmatch{q,4}(1); %N2 region length
%                 CompareMat(q,7) = (calcTDTscore(N1seq)*length(N1seq))^2; %Dmatch{q,4}(3); %N1 region length

                q = q+1;
            end
        end
    
        %Determine maximum alignment score for D
        TotScore = sum(CompareMat(:,3:7),2);%.*CompareMat(:,6).*CompareMat(:,7);
        BestD = TotScore == max(TotScore);    
        if sum(BestD) > 1 %Break ties by looking at the D end deletion counts
            D5D3del = cell2mat(Dmatch(:,3));
            MaxDels = max(D5D3del(:,[1 3]),[],2);
            BestD = BestD & (MaxDels == min(MaxDels(BestD)));
        end    
        BestMatch = find(BestD == 1);
        BestMatch = BestMatch(1); %If still tie, take 1st one only.
        
        if BestMatch == 1; continue; end %no changes needed after all
        %[TotScore(1,end) TotScore(BestMatch,end)]
        %Update the necessary informations for the Tdata
        Dmatch = Dmatch(BestMatch,:);
        VnewDel = CompareMat(BestMatch,1); %Nts to trim from V portion
        JnewDel = CompareMat(BestMatch,2); %Nts to trim from J portion
        VMDNJnew = [VMDNJ(1)-VnewDel  Dmatch{1,4}  VMDNJ(end)-JnewDel];

        if VMDNJnew(2) >0 %Attempt a Vmatch realignment,  using leftover as M's.
            Vnt = [FirstSeq(1:sum(VMDNJnew(1:2))) Vpad];
            AllowedMiss = ceil(0.05 * (length(Vnt) - VpadCt));
            Vmatch = findGeneMatch(Vnt,Vmap,'V',AllowedMiss); %Redo for all V's
            VMDNJnew(2) = Vmatch{4}(3) - length(Vpad);
            VMDNJnew(1) = sum(Vmatch{4}(1:2));
            VnewDel = Vmatch{3}(1,3) - VdelCur; %Nts to trim/add from V
            VmapNum = Vmatch{1};
            Vname = Vmatch{2};
        end

        if VMDNJnew(4)>0 %Attempt a Jmatch realignment, using leftover as N's.
            Jnt = [Jpad FirstSeq(sum(VMDNJnew(1:3))+1:end)];
            AllowedMiss = ceil(MissRate * (length(Jnt) - JpadCt));
            Jmatch = findGeneMatch(Jnt,Jmap,'J',AllowedMiss); %Redo for all J's
            VMDNJnew(5) = sum(Jmatch{4}(2:3));        
            VMDNJnew(4) = length(FirstSeq) - sum(VMDNJnew([1 2 3 5]));
            JnewDel = Jmatch{3}(1,1) - JdelCur; %Nts to trim/add from J
            JmapNum = Jmatch{1};
            Jname = Jmatch{2};
        end
        
        if VMDNJnew(1)*VMDNJnew(3)*VMDNJnew(5) == 0; %Try not to lose the VDJ.
            continue; 
        end 
        if min(VMDNJnew) < 0; %Don't fix errors
            disp(['Weird D correction found in Group ' num2str(GrpNumUnq(y))]);
            continue; 
        end 

        VDDJdels = [(VdelCur+VnewDel)  Dmatch{1,3}(1,1)  Dmatch{1,3}(1,3)  (JdelCur+JnewDel)];
        Tdata(:,DelLoc) = repmat(num2cell(VDDJdels),size(Tdata,1),1);
        Tdata(:,LengthLoc) = repmat(num2cell(VMDNJnew),size(Tdata,1),1);
        Tdata(:,FamNumLoc(2)) = Dmatch(1,1);
        Tdata(:,FamLoc(2)) = Dmatch(1,2);
        Tdata(:,FamNumLoc(1)) = {VmapNum};
        Tdata(:,FamLoc(1)) = {Vname};
        Tdata(:,FamNumLoc(3)) = {JmapNum};
        Tdata(:,FamLoc(3)) = {Jname};

        %Redo the alignments and classifiers
        Tdata = buildVDJalignment(Tdata,NewHeader,Vmap,Dmap,Jmap);
        Tdata = buildRefSeq(Tdata,NewHeader,'same','germline','first');%Same length, germline substitution, on first sequence of each group
        Tdata = makeClassifier(Tdata,NewHeader);
        Tdata = appendMutCt(Tdata,NewHeader); %SHM infor on the VMDNJ segments

        VDJdata(IdxLoc,:) = Tdata;
    end
end