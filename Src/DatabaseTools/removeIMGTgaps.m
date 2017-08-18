%removeIMGTgaps will get a V gene sequence and remove gaps. Will return a
%string data that can be used to get the gaps back.
%
%  [NoGapSeq, GapInfo] = removeIMGTgaps(Seq)
%
%  INPUT
%    Seq: nt sequence given by IMGT
%  
%  OUTPUT
%    NoGapSeq: nt sequence without gaps, replacing wildcard gaps with 'X'
%    GapInfo: A string sequence in the format M(-G(M(-G, where +M numbers
%      are nt counts, and -G numbers are gap counts.
%      count. This is used to return a NoGapSeq to a GappedSeq
%
%  EXAMPLE
%    Seq = 'cag...gtt...act...gaatctgca...ctcac';
%    [NoGapSeq, GapInfo] = removeIMGTgaps(Seq)
%    NoGapSeq =
%      caggttactgaatctgcactcac
%    GapInfo =
%      3#-3#3#-3#3#-3#9#-3#5
%
%  See also addIMGTgaps

function [NoGapSeq, GapInfo] = removeIMGTgaps(Seq)
DotLoc = Seq == '.';
GapTracker = cell(1, length(Seq));
j = 1; %Seq counter
k = 1; %GapTracker counter
while j <= length(Seq)
    if DotLoc(j) == 0
        StartPos = j;
        while DotLoc(j) == 0
            j = j + 1;
            if j > length(Seq); break; end
        end
        GapTracker{k} = j - StartPos; %Postive value for nts
        k = k + 1;
    else
        StartPos = j;
        while DotLoc(j) == 1
            j = j + 1;
            if j > length(Seq); break; end
        end
        GapTracker{k} = -(j - StartPos); %Negative value for gaps
        k = k + 1;
    end
end
GapTracker(k:end) = [];
GapInfo = sprintf(repmat('%d#', 1, length(GapTracker)), GapTracker{:});
GapInfo(end) = [];

NoGapSeq = Seq(~DotLoc);
