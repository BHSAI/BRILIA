%findEquivMatch will open up a VDJdata file, and then find the equivalent
%NT matches possible per VDJ gene based on alignment identity. 100% matches
%are then grouped to yield the multiplicity family name. EX
%if 90 nts of V gene can be exactly matched to 3 V gene family, then ALL 3
%gene become equivalent matches, and thus VDJ family name and number are
%extended to include all 3. 
%
%THis is used mainly for simulated data sets, to get equivalent match
%family data.

function VDJdata = findEquivMatch()
[VDJdata, VDJheader, FileName, FilePath] = openSeqData;
%Open the VDJ Gene database in .mat file, if it exist. 
[Vmap,Dmap,Jmap] = getCurrentDatabase;

H = getHeaderVar(VDJheader);

%Expand the gene family names, using "equivalent" match
%Equivalent if non-mutated nts are the same. For mutants, replace those
%with "X", and then find equivalent matches. 
for j = 1:size(VDJdata,1)
    j
    VMDNJ = cell2mat(VDJdata(j,H.LengthLoc));
    Seq = VDJdata{j,H.SeqLoc};
    SeqV = Seq(1:VMDNJ(1));
    SeqD = Seq(sum(VMDNJ(1:2))+1:sum(VMDNJ(1:3)));
    SeqJ = Seq(end-VMDNJ(end)+1:end);
    
%     SeqRef = VDJdata{j,H.RefSeqLoc};
%     SeqVRef = SeqRef(1:VMDNJ(1));
%     SeqDRef = SeqRef(sum(VMDNJ(1:2))+1:sum(VMDNJ(1:3)));
%     SeqJRef = SeqRef(end-VMDNJ(end)+1:end);
%     
%     Vmiss = SeqVRef ~= SeqV;
%     Dmiss = SeqDRef ~= SeqD;
%     Jmiss = SeqJRef ~= SeqJ;
%     
%     SeqV(Vmiss) = 'X';
%     SeqD(Dmiss) = 'X';
%     SeqJ(Jmiss) = 'X';
%     
%     VcurDel = VDJdata{j,H.DelLoc(1)};
%     D5curDel = VDJdata{j,H.DelLoc(2)};
%     D3curDel = VDJdata{j,H.DelLoc(3)};
%     JcurDel = VDJdata{j,H.DelLoc(4)};

    %======================================================================
    EquivV = zeros(size(Vmap,1),1);
    for k = 1:size(Vmap,1)
        VrefSeq = Vmap{k,1};%(end-VMDNJ(1)-VcurDel+1:end-VcurDel);
        [Score,~, ~] = convolveSeq(SeqV,VrefSeq,0);
        if Score(1) >= 0.98*length(SeqV) || Score(1) >= 0.98*length(VrefSeq);%^2
            EquivV(k) = k;
        end
    end
    EquivV(EquivV==0) = [];
    AllName = '';
    for q = 1:length(EquivV)
        if q == 1
            AllName = Vmap{EquivV(q),3};
        else
            AllName = [AllName '|' Vmap{EquivV(q),3}];    
        end
    end
    VDJdata{j,H.FamNumLoc(1)} =  EquivV;
    VDJdata{j,H.FamLoc(1)} = AllName;
    
    %======================================================================
    EquivD = zeros(size(Dmap,1),1);
    for k = 1:size(Dmap,1)
        DrefSeq = Dmap{k,1};
        if length(DrefSeq) < VMDNJ(3); continue; end
%         if length(DrefSeq) < VMDNJ(3)+D5curDel+D3curDel;
%             continue
%         else
%             DrefSeq = DrefSeq(end-VMDNJ(3)-D3curDel+1:end-D3curDel);
%         end
        [Score,~, ~] = convolveSeq(SeqD,DrefSeq,0);
        if Score(1) >= 0.98*length(SeqD) || Score(1) >= 0.98*length(DrefSeq)%^2
            EquivD(k) = k;
        end
    end
    EquivD(EquivD==0) = [];
    AllName = '';
    for q = 1:length(EquivD)
        if q == 1
            AllName = Dmap{EquivD(q),3};
        else
            AllName = [AllName '|' Dmap{EquivD(q),3}];    
        end
    end
    VDJdata{j,H.FamNumLoc(2)} =  EquivD;
    VDJdata{j,H.FamLoc(2)} = AllName;
    
    %======================================================================
    EquivJ = zeros(size(Jmap,1),1);
    for k = 1:size(Jmap,1)
        JrefSeq = Jmap{k,1};%(1+JcurDel:VMDNJ(end)+JcurDel);
        [Score,~, ~] = convolveSeq(SeqJ,JrefSeq,0);
        if Score(1) >= 0.98*length(SeqJ) || Score(1) >= 0.98*length(JrefSeq)%^2
            EquivJ(k) = k;
        end
    end
    EquivJ(EquivJ==0) = [];
    AllName = '';
    for q = 1:length(EquivJ)
        if q == 1
            AllName = Jmap{EquivJ(q),3};
        else
            AllName = [AllName '|' Jmap{EquivJ(q),3}];    
        end
    end
    VDJdata{j,H.FamNumLoc(3)} =  EquivJ;
    VDJdata{j,H.FamLoc(3)} = AllName;    
end


%Before saving to xlsx, convert columns with matrix values into char for saving
DotLoc = find(FileName == '.');
FileNamePre = FileName(1:DotLoc(end)-1);
for d1 = 1:size(VDJdata,1)
    for d2 = 1:3
        VDJdata{d1,H.FamNumLoc(d2)} = mat2str(VDJdata{d1,H.FamNumLoc(d2)});
    end
end
xlswrite([FileNamePre '_Equiv.xlsx'],[VDJheader; VDJdata]);
