%countConsecLoc will count the largest consecutive 1's or 0's. 
%  EX: MatchResult = '    ||||||  | ||||'
%  countConsecLoc(MatchResults == '|') = 6

function ConsecCount = countConsecLoc(BinData)
ConsecCount = 0;
CurCount = 0;
OneLoc = find(BinData == 1);
for j = 1:length(OneLoc)-1
    if OneLoc(j) + 1 == OneLoc(j+1)
        if CurCount == 0
            CurCount = 2;
        else
            CurCount = CurCount + 1;
        end
    else
        if CurCount > ConsecCount
            ConsecCount = CurCount;
        end
        CurCount = 0;
    end
end
%Perform last iteratino
if CurCount > ConsecCount
    ConsecCount = CurCount;
end
