%calcAromaticity will calculate the num of aromatic residues / total
%residues in a sequence. Aromatic Residues = FWY.
%
%  HPI = calcAromaticity(Seq)
%
%  INPUT
%    Seq: string or cell of string of amino acid sequences
%
%  OUTPUT
%    Aromaticity: value or matrix of Aromaticity
%
%  REFERENCE
%  Lobry, J. R. and C. Gautier (1994). "Hydrophobicity, expressivity and
%  aromaticity are the major trends of amino-acid usage in 999 Escherichia
%  coli chromosome-encoded genes." Nucleic Acids Research 22(15):
%  3174-3180.

function Aromaticity = calcAromaticity(Seq)
Seq = upper(Seq);
if iscell(Seq)
    Aromaticity = cellfun(@(x) countLetter(x, 'FWY'), Seq);
else
    Aromaticity = countLetter(Seq, 'FWY');
end

function Score = countLetter(Seq, Query)
Score = 0;
for j = 1:length(Query)
    Score = Score + sum(Query(j) == Seq);
end
Score = Score / length(Seq);