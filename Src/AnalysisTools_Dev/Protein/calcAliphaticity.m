%calcAliphaticity will calculate the num of aliphatic residues / total
%residues in a sequence. Aliphatic Residues = GAVLI.
%
%  HPI = calcAliphaticity(Seq)
%
%  INPUT
%    Seq: string or cell of string of amino acid sequences
%
%  OUTPUT
%    Aliphaticity: value or matrix of Aliphaticity
%
%  REFERENCE (based on Aromaticity definition here)
%  Lobry, J. R. and C. Gautier (1994). "Hydrophobicity, expressivity and
%  aromaticity are the major trends of amino-acid usage in 999 Escherichia
%  coli chromosome-encoded genes." Nucleic Acids Research 22(15):
%  3174-3180.

function Aliphaticity = calcAliphaticity(Seq)
Seq = upper(Seq);
if iscell(Seq)
    Aliphaticity = cellfun(@(x) countLetter(x, 'GAVLI'), Seq);
else
    Aliphaticity = countLetter(Seq, 'GAVLI');
end

function Score = countLetter(Seq, Query)
Score = 0;
for j = 1:length(Query)
    Score = Score + sum(Query(j) == Seq);
end
Score = Score / length(Seq);
