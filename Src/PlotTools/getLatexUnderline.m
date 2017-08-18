%Simple function for adding underlines to name for latex interpretation.
%
%  LatexName = addLatexUnderline(Str, Position)
%
%  INPUT
%    Str: string
%    Position: the Nth position in Str to underline
%
%  OUTPUT
%    LatexName: string with \underline{} for each position that needs it
%
%  EXAMPLE
%    Str = 'ACGTGTAGTGATG'
%    Position = [1 3:4 6 8]
%    LatexName = getLatexUnderline(Str, Position)
%    LatexName = 
%       \underline{A}C\underline{GT}G\underline{T}A\underline{G}TGATG
function LatexName = getLatexUnderline(Str, Position)
UndIdx = zeros(1, length(Str));
UndIdx(Position) = 1;
j = 1;
c = 1; 
d = 1;
while j <= length(Str)
    if UndIdx(j) > 0
        while UndIdx(j) > 0
            UndIdx(j) = c;
            j = j+1;
            if j > length(Str); break; end
        end
        c = c+1;
    end
    if j > length(Str); break; end
    if UndIdx(j) == 0
        while UndIdx(j) == 0
            UndIdx(j) = -d;
            j = j+1;
            if j > length(Str); break; end
        end
        d = d+1;
    end
end

MaxNum = max(abs(UndIdx));
CellNoUnder = cell(1, MaxNum);
for j = 1:max(abs(UndIdx(UndIdx < 0)))
    CellNoUnder{j} = Str(UndIdx == -j);
end
CellUnder = cell(1, MaxNum);
for j = 1:max(abs(UndIdx(UndIdx > 0)))
    CellUnder{j} = ['\underline{' Str(UndIdx == j) '}'];
end

if UndIdx(1) > 0 %Underscore comes first
    AllCell = [CellUnder; CellNoUnder];
else
    AllCell = [CellNoUnder; CellUnder];
end
AllCell = AllCell(:);
LatexName = sprintf(repmat('%s', 1, length(AllCell)), AllCell{:});
