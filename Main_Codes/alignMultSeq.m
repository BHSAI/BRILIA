%alignMultSeq will take a cell string of sequences, and align them with
%respect to 1st seq in the series
%MinScore = alignment score, 4 consec match = 4^2 = 16. 

function varargout = alignMultSeq(SeqSet,MinScore)
if size(SeqSet,1) == 1;
    varargout{1} = SeqSet;
    if nargout >= 2
        varargout{2} = SeqSet{1};
    end
    return
end

LeftPadCt = zeros(length(SeqSet),2); %LeftPadding Required
LeftPadCt(2:end,2) = 1; %LeftPadding Required
g = 1;
while g <= size(SeqSet,1)
    RefIdx = find(LeftPadCt(:,2) == 0);
    SeqIdx = find(LeftPadCt(:,2) ~= 0);
    BreakThis = 0;
    for r = 1:length(RefIdx)
        RefSeq = SeqSet{RefIdx(r)};
        for s = 1:length(SeqIdx)
            CurSeq = SeqSet{SeqIdx(s)};
            [Score, ~, StartAt] = convolveSeq(RefSeq,CurSeq);
            if Score(3) >= MinScore %consecutive 4 found
                if StartAt(2) < 0
                    LeftPadCt(SeqIdx(s),1) = abs(StartAt(2)) + LeftPadCt(RefIdx(r));
                else
                    LeftPadCt(SeqIdx(s),1) = StartAt(1) - StartAt(2) + LeftPadCt(RefIdx(r));
                end
                LeftPadCt(SeqIdx(s),2) = 0;
                BreakThis = 1;
                break
            end
        end
        if BreakThis == 1
            break
        end        
    end
    g = g+1;
end

%Make sure all negative left pads are converted to minimum 0.
if min(LeftPadCt(:,1)) < 0
    LeftPadCt(:,1) = abs(min(LeftPadCt(:,1))) + LeftPadCt(:,1);
end

%Determine the right size padding, based on maximum length sequence
MaxLen = 0;
for j = 1:size(SeqSet,1)
    CurLen = length(SeqSet{j}) + LeftPadCt(j,1);
    if CurLen > MaxLen
        MaxLen = CurLen;
    end
    LeftPadCt(j,2) = CurLen;
end
LeftPadCt(:,2) = MaxLen - LeftPadCt(:,2);

%Add padding to left and right sides of SeqSet
for j = 1:length(SeqSet)
    SeqSet{j} = [repmat('-',1,LeftPadCt(j,1)) SeqSet{j} repmat('-',1,LeftPadCt(j,2))];
end

%Prepare outputs
varargout{1} = SeqSet;
if nargout >= 2
    CharSeq = char(SeqSet);
    CharSeq(CharSeq == 'B') = 'V'; %Temporarily convert B bonding to Z.
    CharSeq(CharSeq == 'U') = 'Y'; %Temporarily convert B bonding to Z.
    ConsMat = zeros(20,size(CharSeq,2));
    for c = 1:size(CharSeq,2)
        ConsMat(:,c) = cell2mat(struct2cell(aacount(CharSeq(:,c))));
    end
    ConsMat = ConsMat / size(CharSeq,1);
    ConsMax = max(ConsMat,[],1) >= 0.70;
    StartEndLoc = find(ConsMax == 1);
    StartEndLoc = [StartEndLoc(1) StartEndLoc(end)];
    ConsSeq = seqconsensus(CharSeq);
    ConsSeq = ConsSeq(StartEndLoc(1):StartEndLoc(end));
    ConsSeq(ConsSeq == 'V') = 'B';
    ConsSeq(ConsSeq == 'Y') = 'U';
    varargout{2} = ConsSeq;
end
