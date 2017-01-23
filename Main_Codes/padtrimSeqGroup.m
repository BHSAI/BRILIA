%padtrimSeqGroup will trim and pad sequences within a group such that they are
%all the same length. This is required to ensure phylogeny trees are
%constructed correctly and conformGeneGroup works.
%
%  VDJdata = padtrimSeqGroup(VDJdata,NewHeader,GroupBy,KeepMode)
%    INPUT
%      GroupBy ['grpnum','cdr3length']: Will group sequence by group number
%        or CDR3 lengths first, before ensuring each group have same-length
%        sequences.
%      KeepMode ['mean','min','max','trim']: Sets how many nts left of
%        CDR3start and right of CDR3end to keep, based on the mean, min, or
%        max length of nts within each group. 'trim' will trim only 'x' out
%        of sequences edges until a non-x letter is reached.
%
%    NOTE
%      Use GroupBy 'cdr3length' before clustering by phylogeny tree, and
%      then GroupBy 'grpnum' after that is done. Use KeepMode 'max' to
%      retain the full sequences of all groups.
%
%    WARNING
%      This WILL remove entries in VDJdata that have no resolved CDR3start
%      and CDR3end, since these sequences cannot be clustered correctly by
%      same-length CDR3 regions.
%
%  See also clusterGene, conformGeneGroup, padtrimSeq.

function [VDJdata, BadVDJdata] = padtrimSeqGroup(VDJdata,NewHeader,GroupBy,KeepMode,varargin)
getHeaderVar;

%Ensure no empty values in the CDR3Loc column
EmpLoc1 = findCell(VDJdata(:,CDR3Loc(3)),[]);
EmpLoc2 = findCell(VDJdata(:,CDR3Loc(4)),[]);
if min(EmpLoc1) > 0
    VDJdata(EmpLoc1,CDR3Loc(3)) = {0};
end
if min(EmpLoc2) > 0
    VDJdata(EmpLoc2,CDR3Loc(4)) = {0};
end

%Get the length of sequences
SeqLengths = zeros(size(VDJdata,1),1);
for j = 1:length(SeqLengths)
    SeqLengths(j) = length(VDJdata{j,SeqLoc});
end

%Identify start and end locations of CDR3 region
CDR3starts = cell2mat(VDJdata(:,CDR3Loc(3)));
CDR3ends = cell2mat(VDJdata(:,CDR3Loc(4)));

%Reset CDR3ends that go beyond seq length
% ExtendCDR3loc = CDR3ends > SeqLengths;
% CDR3endsT = CDR3ends; %Only for the purpose of getting CDR3length. CDR3ends can be longer than seq length, depending on if it was cut off.
% CDR3endsT(ExtendCDR3loc) = SeqLengths(ExtendCDR3loc);
CDR3Lengths = CDR3ends - CDR3starts + 1;

%Remove those without CDR3s
InvalidLoc = CDR3starts == 0 | CDR3ends == 0 | CDR3Lengths <= 0;
BadVDJdata = VDJdata(InvalidLoc,:);
VDJdata(InvalidLoc,:) = [];
CDR3starts(InvalidLoc) = [];
CDR3ends(InvalidLoc) = [];
CDR3Lengths(InvalidLoc) = [];
if sum(InvalidLoc) > 0
    disp(['Removed ' num2str(sum(InvalidLoc)) ' seqs lacking resolved CDR3s']);
end

%Determine how to group sequences
switch lower(GroupBy)
    case 'grpnum'
        GrpNum = cell2mat(VDJdata(:,GrpNumLoc));
        [~,~,UnqIdx] = unique(GrpNum);
    case 'cdr3length'
        [~,~,UnqIdx] = unique(CDR3Lengths);
end

%Trip/pad sequences by group
UpdateIdx = zeros(size(VDJdata,1),1,'logical');
for y = 1:max(UnqIdx)
    GrpIdx = find(UnqIdx == y);
    
    %Determine how many nts left of C and right of W to keep
    switch KeepMode
        case 'mean'
            KeepLeft = round(mean(CDR3starts(GrpIdx)-1));
            KeepRight = round(mean(SeqLengths(GrpIdx) - CDR3ends(GrpIdx)));
        case 'min'
            KeepLeft = min(CDR3starts(GrpIdx)-1);
            KeepRight = min(SeqLengths(GrpIdx) - CDR3ends(GrpIdx));
        case 'max'
            KeepLeft = max(CDR3starts(GrpIdx)-1);
            KeepRight = max(SeqLengths(GrpIdx) - CDR3ends(GrpIdx));
        case 'trim'
            %Look for longest distance from CDR3point to non-x int
            DelLeft = inf;
            DelRight = inf;
            for j = 1:length(GrpIdx)
                Seq = VDJdata{GrpIdx(j),SeqLoc};
                NonXloc = regexpi(Seq,'[^X*]');
                CurDelLeft = NonXloc(1)-1;
                CurDelRight = length(Seq)-NonXloc(end);
                if CurDelLeft < DelLeft
                    DelLeft = CurDelLeft;
                end
                if CurDelRight < DelRight
                    DelRight = CurDelRight;
                end
            end
    end
    
    %Perform trim/pad
    for j = 1:length(GrpIdx)
        Seq = VDJdata{GrpIdx(j),SeqLoc};
        RefSeq = VDJdata{GrpIdx(j),RefSeqLoc};
        Vlen = VDJdata{GrpIdx(j),LengthLoc(1)};
        Jlen = VDJdata{GrpIdx(j),LengthLoc(5)};
        CDR3start = CDR3starts(GrpIdx(j));
        CDR3end = CDR3ends(GrpIdx(j));
        
        if strcmpi(KeepMode,'trim') %Trim only, special function here.
            %Update seq and variables of importance
            NewSeqNT = Seq(DelLeft+1:end-DelRight);
            NewRefSeqNT = RefSeq(DelLeft+1:end-DelRight);
            NewVlen = Vlen - DelLeft;
            NewJlen = Jlen - DelRight;
            NewCDR3start = CDR3start - DelLeft;
            NewCDR3end = CDR3end - DelLeft;
        else
        
            %Determine the new start loc of seq
            StartLoc = CDR3starts(GrpIdx(j)) - KeepLeft;
            if StartLoc <= 0 %Need padding X
                LeftPad = repmat('X',1,1 - StartLoc);
                StartLoc = 1;
                CDR3shift = length(LeftPad);
            else
                LeftPad = '';
                CDR3shift = -StartLoc + 1;
            end

            %Determine the new end loc of seq
            EndLoc = CDR3ends(GrpIdx(j)) + KeepRight;
            if EndLoc > length(Seq) %need padding X
                RightPad = repmat('X',1,EndLoc - length(Seq));
                EndLoc = length(Seq);
            else
                RightPad = '';
            end

            %Update seq and variables of importance
            NewSeqNT = [LeftPad, Seq(StartLoc:EndLoc), RightPad];
            NewRefSeqNT = [LeftPad, RefSeq(StartLoc:EndLoc), RightPad];
            NewVlen = Vlen + length(LeftPad) - StartLoc + 1;
            NewJlen = Jlen + length(RightPad) - (length(Seq) - EndLoc);
            NewCDR3start = CDR3start + CDR3shift;
            NewCDR3end = CDR3end + CDR3shift;
        end
        
        VDJdata(GrpIdx(j),[SeqLoc RefSeqLoc LengthLoc(1) LengthLoc(5) CDR3Loc(3) CDR3Loc(4)]) = {NewSeqNT NewRefSeqNT NewVlen NewJlen NewCDR3start NewCDR3end};
        UpdateIdx(GrpIdx(j)) = 1;
    end
end

%Update those that have changed
VDJdata(UpdateIdx,:) = updateVDJdata(VDJdata(UpdateIdx,:),NewHeader,varargin);
