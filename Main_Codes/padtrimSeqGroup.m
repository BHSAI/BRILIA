%padtrimSeqGroup will trim and pad sequences within a group such that they
%are all the same length. This is required to ensure phylogeny trees are
%constructed correctly and conformGeneGroup works.
%
%  VDJdata = padtrimSeqGroup(VDJdata,NewHeader,GroupBy,KeepMode,BasedOn);
%
%  VDJdata = padtrimSeqGroup(...,Vmap,Dmap,Jmap);
%
%    INPUT
%      GroupBy ['grpnum','cdr3length']: Will group sequence by group number
%        or CDR3 lengths first, before ensuring each group have same-length
%        sequences.
%      KeepMode ['mean','min','max','trim']: Sets how many nts left of
%        CDR3start and right of CDR3end to keep, based on the mean, min, or
%        max length of nts within each group. 'trim' will trim only 'x' out
%        of sequences edges until a non-x letter is reached.
%      BasedOn ['Seq','RefSeq']: determines how to trim sequences based on
%        either the flanking X nts in Seq or RefSeq. Use 'RefSeq' if you
%        have sequences that are longer than V and J germline genes.
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

function [VDJdata, BadVDJdata] = padtrimSeqGroup(VDJdata,NewHeader,GroupBy,KeepMode,BasedOn,varargin)
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
CDR3Lengths = CDR3ends - CDR3starts + 1;

%Remove those without CDR3s
InvalidLoc = (CDR3starts <= 0) | (CDR3ends <= 0) | (CDR3Lengths <= 3) | (CDR3ends > SeqLengths) | (CDR3starts > SeqLengths); %Ex, length of 3, being TGT, cannot be both CDR3start to CDR3end.
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

%Determine which one to base trimming on
if strcmpi(KeepMode,'trim') && strcmpi(BasedOn,'RefSeq');
    EvalSeqLoc = RefSeqLoc; %Uses flanking X's in RefSeq to decide trim
    OtherSeqLoc = SeqLoc; %Will trim Seq to keep the lengths the same
else
    EvalSeqLoc = SeqLoc; %Uses flanking X's in Seq ot decide trim
    OtherSeqLoc = RefSeqLoc; %Will trim RefSeq to keep the lenghts the same
end

%Trim/pad sequences of a group
UpdateIdx = zeros(size(VDJdata,1),1,'logical');
for y = 1:max(UnqIdx)
    GrpIdx = find(UnqIdx == y);    
    CDR3s = CDR3starts(GrpIdx);
    SeqLen = SeqLengths(GrpIdx);

    %Determine how many nts left and right of 104C's 1st nt to keep.
    MinDelLeft = 0;
    MinDelRight = 0;
    switch KeepMode
        case 'mean'
            KeepLeft = round(mean(CDR3s - 1));
            KeepRight = round(mean(SeqLen - CDR3s));
        case 'min'
            KeepLeft = min(CDR3s - 1);
            KeepRight = min(SeqLen - CDR3s);
        case 'max'
            KeepLeft = max(CDR3s - 1);
            KeepRight = max(SeqLen - CDR3s);
        case 'trim'
            KeepLeft = max(CDR3s - 1);
            KeepRight = max(SeqLen - CDR3s);

            %Fill up the left and right X deletions
            MinDelLeft = max(SeqLen);
            MinDelRight = max(SeqLen);
            for j = 1:length(GrpIdx)
                EvalSeq = VDJdata{GrpIdx(j),EvalSeqLoc};
                NonXloc = regexpi(EvalSeq,'[^X*]');
                CurDelLeft = NonXloc(1)-1;
                CurDelRight = length(EvalSeq)-NonXloc(end);
                if CurDelLeft < MinDelLeft
                    MinDelLeft = CurDelLeft;
                end
                if CurDelRight < MinDelRight
                    MinDelRight = CurDelRight;
                end
                
                %Do first one only, if using RefSeq since that's germline
                if strcmpi(BasedOn,'RefSeq')
                    break
                end
            end
    end
    
    %Calculate how much left / right to pad(+ value) or trim(- value)
    LeftAdd = KeepLeft - (CDR3s - 1) - MinDelLeft;
    RightAdd = KeepRight + CDR3s - SeqLen - MinDelRight;
    
    %Perform trim/pad
    for j = 1:length(GrpIdx)
        %Don't update if there's nothing to update
        if (LeftAdd(j) == 0) && (RightAdd(j) == 0)
            continue
        end

        %Extract variables to update
        EvalSeq = VDJdata{GrpIdx(j),EvalSeqLoc};
        OtherSeq = VDJdata{GrpIdx(j),OtherSeqLoc};
        Vlen = VDJdata{GrpIdx(j),LengthLoc(1)};
        Jlen = VDJdata{GrpIdx(j),LengthLoc(5)};
        CDR3start = CDR3starts(GrpIdx(j));
        CDR3end = CDR3ends(GrpIdx(j));
        
        %Determine left side edits
        if LeftAdd(j) > 0
            LeftPad = repmat('X',1,LeftAdd(j));
            S1 = 1;
        else
            LeftPad = '';
            S1 = abs(LeftAdd(j)) + 1;
        end
        
        %Determine right side edits
        if RightAdd(j) > 0
            RightPad = repmat('X',1,RightAdd(j));
            S2 = length(EvalSeq);
        else
            RightPad = '';
            S2 = length(EvalSeq) + RightAdd(j); %Remember RightAdd is neg for trim
        end
        
        %Update seq and variables of importance
        NewEvalSeq = [LeftPad, EvalSeq(S1:S2), RightPad];
        NewOtherSeq = [LeftPad, OtherSeq(S1:S2), RightPad];
        NewVlen = Vlen + LeftAdd(j);
        NewJlen = Jlen + RightAdd(j);
        NewCDR3start = CDR3start + LeftAdd(j);
        NewCDR3end = CDR3end + LeftAdd(j);
        
        VDJdata(GrpIdx(j),[EvalSeqLoc OtherSeqLoc LengthLoc(1) LengthLoc(5) CDR3Loc(3) CDR3Loc(4)]) = {NewEvalSeq NewOtherSeq NewVlen NewJlen NewCDR3start NewCDR3end};
        UpdateIdx(GrpIdx(j)) = 1;
    end
end

%Update those that have changed
VDJdata(UpdateIdx,:) = updateVDJdata(VDJdata(UpdateIdx,:),NewHeader,varargin);
