%makeClassifier will make a letter sequence that labels each nt by the
%segment it belongs to. The "Classifier" sequence is stored in VDJdata in
%the "Classifier" column.
%
%  VDJdata = makeClassifier(VDJdata,NewHeader)
%
%  EXAMPLE
%    Classifier = 'VVVVVVVVVVVVVVVVppmmmppDDDDDDppnnppJJJJJJJJJ';
%      V = V   segment
%      M = Nvd region
%      D = D   segment
%      N = Ndj region
%      J = J   segment
%      B = p-nucleotides in the Nvd region
%      P = p-nucleotides in the Ndj region
%      v,m,d,n,j,b,p = lowercase are mismatched with the reference sequence

function VDJdata = makeClassifier(VDJdata,NewHeader)
getHeaderVar;

%Correct the classifiers according to the group, mainly p and b's only.
GrpNum = cell2mat(VDJdata(:,GrpNumLoc));
UnqGrpNum = unique(GrpNum);
for y = 1:length(UnqGrpNum)
    try
        IdxLoc =  find(UnqGrpNum(y) == GrpNum);
        j = IdxLoc(1);

        %Extract information
        Seq = VDJdata{j,SeqLoc};
        RefSeq = VDJdata{j,RefSeqLoc};
        VMDNJ = cell2mat(VDJdata(j,LengthLoc));

        if sum(VMDNJ) ~= length(Seq) || sum(VMDNJ) ~= length(RefSeq)
            continue;
        end

        %Assemble the starting point classifier
        ClassV = repmat('V',1,VMDNJ(1)); %V segment
        ClassM = repmat('M',1,VMDNJ(2)); %Nvd segment
        ClassD = repmat('D',1,VMDNJ(3)); %D segment
        ClassN = repmat('N',1,VMDNJ(4)); %Ndj segment
        ClassJ = repmat('J',1,VMDNJ(5)); %J segment

        %To find p-nucleotides, see if germline deletions are 0.
        VrefDel3 = VDJdata{j,DelLoc(1)};
        DrefDel5 = VDJdata{j,DelLoc(2)};
        DrefDel3 = VDJdata{j,DelLoc(3)};
        JrefDel5 = VDJdata{j,DelLoc(4)};

        RefClass = sprintf('%s%s%s%s%s',ClassV,ClassM,ClassD,ClassN,ClassJ);

        %Find the Nvd and Ndj p-nts.
        P_V3 = findPnts(RefSeq,VMDNJ(1),'right',0,VrefDel3);
        P_D5 = findPnts(RefSeq,sum(VMDNJ(1:2))+1,'left',DrefDel5,0); 
        P_D3 = findPnts(RefSeq,sum(VMDNJ(1:3)),'right',0,DrefDel3);
        P_J5 = findPnts(RefSeq,sum(VMDNJ(1:4))+1,'left',JrefDel5,0);
        
        %Label Nvd p-nts as B, and Ndj p-nts as P
        RefClass(P_V3 | P_D5) = 'B';
        RefClass(P_D3 | P_J5) = 'P';

        %Fillin the Classifier and Formatted Seq 
        for k = 1:length(IdxLoc)
            CurSeq = VDJdata{IdxLoc(k),SeqLoc};
            MissLoc = CurSeq ~= RefSeq;
            SeqClass = RefClass;
            SeqClass(MissLoc) = lower(SeqClass(MissLoc));

            VDJdata{IdxLoc(k),FormClassLoc(1)} = formatSeq(VDJdata{IdxLoc(k),SeqLoc},SeqClass);
            VDJdata{IdxLoc(k),FormClassLoc(2)} = SeqClass;
        end
    catch
        ErrorMsg = sprintf('Errored at %s, sequence # %d',mfilename,y);
        disp(ErrorMsg);
        VDJdata(IdxLoc,MiscLoc) = repmat({ErrorMsg},length(IdxLoc),1);
    end
    clear ClassV ClassM ClassD ClassN ClassJ
end