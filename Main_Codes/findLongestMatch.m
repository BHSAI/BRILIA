%findLongestMatch finds the longest match in an aligned sequence. 
%
%  [Longest Loc] = findLongestMatch(AlignedSeq) takes an aligned sequence
%  from swalign or nwalign and then looks for longest concecutive matches,
%  which includes : and |. Output is the # of consec match "Longest" and
%  the location start and end of that segment "Loc". If there are ties, it
%  will return the last occurrance.

function [LongestConn, Loc] = findLongestMatch(AlignedSeq)
%Checking input for location of the alignment markers.
if ischar(AlignedSeq)
    if size(AlignedSeq,1) == 3;
        ConnComp = AlignedSeq(2,:) ~= ' ';
    elseif size(AlignedSeq,1) == 1;
        ConnComp = AlignedSeq(1,:) ~= ' ';
    else
        error('Wrong dimension. Must be a char alignement with "|   |||  |"');
    end
elseif isnumeric(AlignedSeq) && min(size(AlignedSeq)) == 1
    ConnComp = AlignedSeq;
elseif islogical(AlignedSeq) && min(size(AlignedSeq)) == 1
    ConnComp = AlignedSeq;
end

%Search for the location and size of longest consecutive match.
LongestConn = 0;
TempLongest = 0;
TempLoc = 1;
for kk = 1:length(ConnComp)
    if ConnComp(kk)
        TempLongest = TempLongest + 1;
        if TempLongest > LongestConn
            LongestConn = TempLongest;
            TempLoc = kk;
        end
    else
        TempLongest = 0;
    end
end

Loc = TempLoc - LongestConn + 1;
