%countCDR3LenFreq will count the CDR3 lengths per clonotype.
%
%  CDR3LenFreq =countCDR3LenFreq
%
%  INPUT
%    VDJdata: BRILIA data cell
%    VDJheader: BRILIA header cell
%
%  OUTPUT
%    CDR3LenFreq: Mx2 cell where 
%      Col 1 - cell of categories as 'xCDR3-N' (N = CDR3 length, x = chain
%              h or l
%      Col 2 - frequency count per category

function CDR3LenFreq = countCDR3LenFreq(VDJdata, VDJheader)
G = getGrpIdx(VDJdata, VDJheader);
Idx = arrayfun(@(x) x.Idx(1), G); 
Map = getVDJmapper(VDJheader);

TempFreq = cell(1, length(Map.Chain));
for c = 1:length(Map.Chain)
    CDR3Len = cell2mat(VDJdata(Idx, Map.([lower(Map.Chain(c)) 'CDR3'])(2)));
    UnqCDR3Len = unique(CDR3Len);
    Freq = cell(length(UnqCDR3Len), 2);
    for k = 1:length(Freq)
        Freq{k, 1} = sprintf('%sCDR-%d', lower(Map.Chain(c)), UnqCDR3Len(k));
        Freq{k, 2} = sum(UnqCDR3Len(k) == CDR3Len);
    end 
    TempFreq{c} = Freq;
end
CDR3LenFreq = vertcat(TempFreq{:});