%findVvsCDR3mut will look at each sequence, find the annotated V framework
%mutation % vs CDR3 region mutations %. Plot of Vmut vs CDR3mut will reveal
%differences in predicted SHM across different annotation results.
%
%Output = [VframeMutFr CDR3MutFr]
%
%  VframeMutFr = nts including the 104C codon / length of VCDR3-Framework
%  length
%
%  CDR3MutFr = nts from 105-117 CDR3 / length of CDR3 regions

function Output = findVvsCDR3mut(varargin)
if isempty(varargin)
    [VDJdata,VDJheader,FileName,FilePath] = openSeqData;
else
    VDJdata = varargin{1};
    VDJheader = varargin{2};
end
H = getHeaderVar(VDJheader);

GrpNum = cell2mat(VDJdata(:,H.GrpNumLoc));
UnqGrpNum = unique(GrpNum);
Output = zeros(size(VDJdata,1),2);
for y = 1:length(UnqGrpNum)
    IdxLoc = find(UnqGrpNum(y) == GrpNum);
    RefSeq = VDJdata{IdxLoc(1),H.SeqLoc}; %Always with respect to root seq.
    for j = 1:length(IdxLoc)
        Seq = VDJdata{IdxLoc(j),H.SeqLoc};
        CDR3loc = cell2mat(VDJdata(IdxLoc(j),H.CDR3Loc(3:4))) + [+3 -3];%Excludes the 104C and 118W
        if CDR3loc(2) > length(Seq);
            CDR3loc(2) = length(Seq);
        end
        VseqFrame = Seq(1:CDR3loc(1)-1);
        VrefFrame = RefSeq(1:CDR3loc(1)-1);
        CDR3seq = Seq(CDR3loc(1):CDR3loc(2));
        CDR3ref = RefSeq(CDR3loc(1):CDR3loc(2));
        Output(IdxLoc(j),1) = sum(VseqFrame ~= VrefFrame)/length(VseqFrame);
        Output(IdxLoc(j),2) = sum(CDR3seq ~= CDR3ref)/length(CDR3seq);
    end
end
% end
% 
% 
% 
% 
% 
% 
% [Vmap,~,Jmap] = getCurrentDatabase;
% 
% %Extract the maximum mutations per clonal group in V_before104C and CDR3
% GrpNum = cell2mat(VDJdata(:,H.GrpNumLoc));
% UnqGrpNum = unique(GrpNum);
% Vmut = zeros(length(UnqGrpNum),1);
% CDR3mut = zeros(length(UnqGrpNum),1);
% for y = 1:length(UnqGrpNum)
%     GrpIdx = find(UnqGrpNum(y) == GrpNum);
%     Tdata = VDJdata(GrpIdx,:);    
%     RefSeq = Tdata{1,H.RefSeqLoc}; %Assumes first seq in group is root
%     
%     %In case there's an issue with RefSeq determination
%     if length(RefSeq) ~= length(VDJdata{GrpIdx(1),H.SeqLoc})
%         continue
%     end
%     
%     %Try to find maximum CDR3 hamming distance between root and seq
%     VmutT = zeros(size(Tdata,1),1);
%     CDR3mutT = zeros(size(Tdata,1),1);
%     for j = 1:size(Tdata,1)
%         SamSeq = Tdata{j,H.SeqLoc};
%         VMDNJ = cell2mat(Tdata(j,H.LengthLoc));
% 
%         %Identify 104C location
%         Vnum = Tdata{j,H.FamNumLoc(1)}(1);
%         VlocC = Vmap{Vnum,end};
%         Vdel = Tdata{j,H.DelLoc(1)};
%         Loc104C = VMDNJ(1) + Vdel - VlocC + 1;
% 
%         %Identify 118W location
%         Jnum = Tdata{j,H.FamNumLoc(3)}(1);
%         JlocW = Jmap{Jnum,end};
%         Jdel = Tdata{j,H.DelLoc(end)};
%         Loc118W = sum(VMDNJ(1:4)) - Jdel + JlocW + 2;
%         if Loc118W > sum(VMDNJ)
%             Loc118W = sum(VMDNJ);
%         end
%         
%         %Compare the sequences
%         VrefFrame = RefSeq(1:Loc104C-1);
%         CDR3refSeq = RefSeq(Loc104C:Loc118W);
% 
%         VsamFrame = SamSeq(1:Loc104C-1);
%         CDR3samSeq = SamSeq(Loc104C:Loc118W);
%         
%         %From nt 1 to 104C before location
%         VmutT(j) = sum(VsamFrame~=VrefFrame)/length(VrefFrame);
%         %From nt 104C start to end of CDR3 118W loc
%         CDR3mutT(j) = sum(CDR3refSeq~=CDR3samSeq)/length(CDR3refSeq);    
%     end
%     
%     %Find when CDR3mutT is max
%     MaxLoc = find(CDR3mutT == max(CDR3mutT));
%     Vmut(y) = VmutT(MaxLoc(1));
%     CDR3mut(y) = CDR3mutT(MaxLoc(1));
% end
% 
% Output = [UnqGrpNum Vmut CDR3mut];
