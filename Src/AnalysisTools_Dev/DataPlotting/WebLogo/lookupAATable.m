%lookupAATable will return the amino acid with the specified HPY and VDWV
%volume.
%
%  AA = loookupAATable(HPI, VDWV)
%
%  INPUT
%    HPI: Hydropathy index to search
%    VDWV: Vad der Waal Volume in A^3 to search
%  
%  OUTPUT
%    AA: amino acid residue that matches exactly to HPI and VDWV

function AminoAcid = lookupAATable(HPI,VDWV)
AATable = getAATable;
AA =  cell2mat(AATable(:,2));
AA_HPY = cell2mat(AATable(:,3));
AA_VDWV = cell2mat(AATable(:,4));
Check1 = AA_HPY == HPI;
Check2 = AA_VDWV == VDWV;
AminoAcid = AA(Check1&Check2);