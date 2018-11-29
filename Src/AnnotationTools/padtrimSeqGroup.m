%padtrimSeqGroup will trim and pad sequences within a group such that they
%are all the same length. This is required to ensure phylogeny trees are
%constructed correctly and conformGeneGroup works.
%
%  [VDJdata, BadVDJdata] = padtrimSeqGroup(VDJdata, Map, GroupBy, KeepMode, BasedOn);
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    Map: Map of the BRILIA data (getVDJmapper(VDJheader))
%    GroupBy ['grpnum', 'cdr3length']: Will group sequence by group number
%      or CDR3 lengths first, before ensuring each group have same-length
%      sequences.
%    KeepMode ['mean', 'min', 'max', 'trim']: Sets how many nts left of
%      CDR3start and right of CDR3end to keep, based on the mean, min, or
%      max length of nts within each group. 'trim' will trim only 'x' out
%      of sequences edges until a non-x letter is reached.
%    BasedOn ['Seq', 'RefSeq']: determines how to trim sequences based on
%      either the flanking X nts in Seq or RefSeq. Use 'RefSeq' if you
%      have sequences that are longer than V and J germline genes.
%
%  OUPUT
%    VDJdata: modified VDJdata with padded/trimmed Seq and RefSeq such that
%      all groups have same-length sequences. Padding/Trimming decisions
%      are "based on" either the RefSeq or Seq. Padding uses "X".
%    BadVDJdata: bad entries of VDJdata that were removed due to lack of
%      resolved CDR3 or info for padding/trimming.
%
%  NOTE
%    Use GroupBy 'cdr3length' before clustering by phylogeny tree, and
%    then GroupBy 'grpnum' after that is done. Use KeepMode 'max' to
%    retain the full sequences of all groups.
%
%  WARNING
%    This WILL remove entries in VDJdata that have no resolved CDR3start
%    and CDR3end, since these sequences cannot be clustered correctly by
%    same-length CDR3 regions.

function [VDJdata, BadVDJdata] = padtrimSeqGroup(VDJdata, Map, GroupBy, KeepMode, BasedOn)
BadVDJdata = {};
if size(VDJdata, 1) == 0; return; end

%Pad or trim Seq and RefSeq based on one of the two sequence type lengths
UpdateIdx = zeros(size(VDJdata, 1), 1, 'logical');
for k = 1:length(Map.Chain)
    C = lower(Map.Chain(k));
    GrpNumLoc = Map.GrpNum;
    SeqLoc    = Map.([C 'Seq']);
    RefSeqLoc = Map.([C 'RefSeq']);
    LengthLoc = Map.([C 'Length']);
    CDR3Loc   = Map.([C 'CDR3']);
    
    %Replace any empty or 1xN CDR3s/e column with 0
    for j = 1:size(VDJdata, 1)
        for b = 3:4
            if isempty(VDJdata{j, CDR3Loc(b)}) || length(VDJdata{j, CDR3Loc(b)}) > 1
                VDJdata{j, CDR3Loc(b)} = 0;
            end
        end
    end

    %Identify beginning and end locations of CDR3
    SeqLens = cellfun('size', VDJdata(:, SeqLoc), 2);
    CDR3Bgns = cell2mat(VDJdata(:, CDR3Loc(3))); %NT pos
    CDR3Ends = cell2mat(VDJdata(:, CDR3Loc(4))); %NT pos
    CDR3Lens = (CDR3Ends - CDR3Bgns + 1)/3;      %AA len

    %Remove invalid CDR3s
    InvalidLoc = (CDR3Bgns <= 1) | (CDR3Ends <= 1) | (CDR3Lens <= 1) | (CDR3Ends > SeqLens) | (CDR3Bgns >= CDR3Ends);
    if any(InvalidLoc)
        BadVDJdata = cat(1, BadVDJdata, VDJdata(InvalidLoc, :));
        VDJdata(InvalidLoc, :) = [];
        CDR3Bgns(InvalidLoc) = [];
        CDR3Ends(InvalidLoc) = [];
        CDR3Lens(InvalidLoc) = [];
        SeqLens(InvalidLoc)  = [];
        fprintf('Removed %d sequences lacking resolved CDR3\n', sum(InvalidLoc));
        if size(VDJdata, 1) == 0; continue; end
    end

    %Determine how to group sequences
    if startsWith(GroupBy, 'grpnum', 'ignorecase', true)
        [~, ~, UnqIdx] = unique(cell2mat(VDJdata(:, GrpNumLoc)));
    elseif startsWith(GroupBy, 'cdr3len', 'ignorecase', true)
        [~, ~, UnqIdx] = unique(CDR3Lens);
    else
        error('%s: Unknown GroupBy option "%s"', mfilename, GroupBy);
    end

    %Determine which seq to base trimming on, Seq or RefSeq
    if strcmpi(KeepMode, 'trim') && strcmpi(BasedOn, 'RefSeq')
        EvalSeqLoc = RefSeqLoc;  %Uses flanking X's in RefSeq to decide trim
        OtherSeqLoc = SeqLoc;    %Will trim Seq to keep the lengths the same
    else
        EvalSeqLoc = SeqLoc;     %Uses flanking X's in Seq to decide trim
        OtherSeqLoc = RefSeqLoc; %Will trim RefSeq to keep the lengths the same
    end

    %Trim/pad sequences of a group
    for y = 1:max(UnqIdx)
        GrpIdx = find(UnqIdx == y);    
        CDR3Bgn = CDR3Bgns(GrpIdx); %This is used for anchoring
        SeqLen = SeqLens(GrpIdx);
        
        %Determine how many nts left and right of 104C's 1st nt to keep.
        MinDelLeft = 0;
        MinDelRight = 0;
        switch KeepMode
            case 'mean'
                KeepLeft = round(mean(CDR3Bgn - 1));
                KeepRight = round(mean(SeqLen - CDR3Bgn));
            case 'min'
                KeepLeft = min(CDR3Bgn - 1);
                KeepRight = min(SeqLen - CDR3Bgn);
            case 'max'
                KeepLeft = max(CDR3Bgn - 1);
                KeepRight = max(SeqLen - CDR3Bgn);
            case 'trim'
                KeepLeft = max(CDR3Bgn - 1);
                KeepRight = max(SeqLen - CDR3Bgn);

                %Fill up the left and right X deletions
                MinDelLeft = max(SeqLen);
                MinDelRight = max(SeqLen);
                for j = 1:length(GrpIdx)
                    EvalSeq = VDJdata{GrpIdx(j), EvalSeqLoc};
                    NonXloc = regexpi(EvalSeq, '[^N*]');
                    CurDelLeft = NonXloc(1) - 1;
                    CurDelRight = length(EvalSeq) - NonXloc(end);
                    if CurDelLeft < MinDelLeft
                        MinDelLeft = CurDelLeft;
                    end
                    if CurDelRight < MinDelRight
                        MinDelRight = CurDelRight;
                    end

                    %Do first one only, if using RefSeq since that's germline
                    if strcmpi(BasedOn, 'RefSeq')
                        break
                    end
                end
        end

        %Calculate how much left / right to pad(+ value) or trim(- value)
        LeftAdd = KeepLeft - (CDR3Bgn - 1) - MinDelLeft;
        RightAdd = KeepRight + CDR3Bgn - SeqLen - MinDelRight;

        %Perform trim/pad
        for j = 1:length(GrpIdx)
            %Don't update if there's nothing to update
            if (LeftAdd(j) == 0) && (RightAdd(j) == 0)
                continue
            end

            %Extract variables to update
            EvalSeq = VDJdata{GrpIdx(j), EvalSeqLoc};
            OtherSeq = VDJdata{GrpIdx(j), OtherSeqLoc};
            Vlen = VDJdata{GrpIdx(j), LengthLoc(1)};
            Jlen = VDJdata{GrpIdx(j), LengthLoc(end)};
            CDR3B = CDR3Bgns(GrpIdx(j));
            CDR3E = CDR3Ends(GrpIdx(j));

            %Determine left side edits
            if LeftAdd(j) >= 0
                LeftPad = repmat('N', 1, LeftAdd(j));
                S1 = 1;
            else
                LeftPad = '';
                S1 = abs(LeftAdd(j)) + 1;
            end

            %Determine right side edits
            if RightAdd(j) >= 0
                RightPad = repmat('N', 1, RightAdd(j));
                S2 = length(EvalSeq);
            else
                RightPad = '';
                S2 = length(EvalSeq) + RightAdd(j); %Remember RightAdd is neg for trim
            end

            %Update seq and variables of importance
            try
                NewEvalSeq  = [LeftPad, EvalSeq(S1:S2),  RightPad];
                NewOtherSeq = [LeftPad, OtherSeq(S1:S2), RightPad];
            catch
                fprintf('%s: ERROR:\n   LeftPad = "%s", RightPad = "%s"\n   S1 = %d, S2 = %d\n   EvalSeq = "%s"\n   OtherSeq = "%s"\n', mfilename, LeftPad, RightPad, S1, S2, EvalSeq, OtherSeq);
                continue
            end
            NewVlen = Vlen + LeftAdd(j);
            NewJlen = Jlen + RightAdd(j);
            NewCDR3B = CDR3B + LeftAdd(j);
            NewCDR3E = CDR3E + LeftAdd(j);

            VDJdata(GrpIdx(j), [EvalSeqLoc OtherSeqLoc LengthLoc(1) LengthLoc(end) CDR3Loc(3) CDR3Loc(4)]) = {NewEvalSeq NewOtherSeq NewVlen NewJlen NewCDR3B NewCDR3E};
            UpdateIdx(GrpIdx(j)) = 1;
        end
    end
end

%Update only SHMs, as all other updates are handled already
VDJdata(UpdateIdx, :) = countSHM(VDJdata(UpdateIdx, :), Map);