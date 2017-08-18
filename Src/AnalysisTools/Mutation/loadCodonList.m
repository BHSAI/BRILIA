%loadCodonList will return a 64x1 cell that lists all 64 codons
%
%  CodonList = loadCodonList
%
%  [CodonList, IntCodonList] = loadCodonList
%
%  OUTPUT
%    CodonList: 64x1 codon list
%    IntCodonList: 64x1 integer codon list (1 = A, 2 = C, 3 = G, 4 = T)
%
function [CodonList, varargout] = loadCodonList()
%Generate 64 codons and it's silent/replacement mutation matrix SRM 
IntCodonList = cell(64, 1);
CodonList = cell(64, 1);
j = 1;
for x = 1:4
    for y = 1:4
        for z = 1:4
            IntCodonList{j} = [x y z];
            CodonList{j} = int2nt([x y z]);
            j = j + 1;
        end
    end
end

if nargout == 2
    varargout{1} = IntCodonList;
end
