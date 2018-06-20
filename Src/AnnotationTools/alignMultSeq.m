%alignMultSeq will take a cell string of sequences, and align them with
%respect to 1st seq in the series
%
%  AlignSeq = alignMultSeq(SeqSet)
%
%  [AlignSeq, ConsSeq] = alignMultSeq(SeqSet)
%
%  INPUT
%    SeqSet: a Mx1 cell of sequences to be aligned, where the 1st sequence
%      is the reference sequence
%
%  OUTPUT
%    AlignSeq: Mx1 cell of sequences that were aligned to the 1st sequence
%      of SeqSet. Paddings are added as '-' to ensure sequences have the
%      same lengths.

function varargout = alignMultSeq(SeqSet)
%If SeqSet is a single sequence, then just return itself
if size(SeqSet,1) == 1
    varargout{1} = SeqSet;
    if nargout >= 2
        varargout{2} = SeqSet{1};
    end
    return
end

%Determine the seq length first
LengthCt = zeros(length(SeqSet),1);
for j = 1:length(SeqSet)
    LengthCt(j) = length(SeqSet{j});
end

%Find left padding requirement
LeftPadCt = zeros(length(SeqSet),2);
SeqA = SeqSet{1};
for b = 2:length(SeqSet)
    SeqB = SeqSet{b};
    [~, StartAt] = alignSeqMEX(SeqA,SeqB, 0, 'r');
    if StartAt(2) < 0 %Need to pad SeqA
        LeftPadCt(b,1) = abs(StartAt(2));
    elseif StartAt(2) > 1 %need to pad SeqB
        LeftPadCt(b,2) = StartAt(2) - 1;
    end
end
LeftPadCt(1,2) = max(LeftPadCt(:,1));
LeftPadCt(2:end,2) = (LeftPadCt(2:end,2) + LeftPadCt(1,2)) - LeftPadCt(2:end,1);
LeftPadCt = LeftPadCt(:,2);

%Find right padding requirement
RightPadCt = max(LengthCt + LeftPadCt) - (LengthCt + LeftPadCt);

%Add padding to left and right sides of SeqSet
for j = 1:length(SeqSet)
    SeqSet{j} = [repmat('-',1,LeftPadCt(j)) SeqSet{j} repmat('-',1,RightPadCt(j))];
end

%Prepare outputs
if nargout >= 1
    varargout{1} = SeqSet;
    if nargout >= 2
        ConsSeq = repmat('-',1,length(SeqSet{1}));
        CharSeq = char(SeqSet);
        for j = 1:length(SeqSet{1})
            TempSeq = strrep(CharSeq(:,j)','-','');
            ConsSeq(j) = mode(TempSeq);
        end
        varargout{2} = ConsSeq;
    end
end
