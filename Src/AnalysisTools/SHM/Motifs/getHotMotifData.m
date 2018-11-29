%getHotMotifData(MotifData) will extract motif and the corresponding
%mutation frequencies for known  hotspots. Returns a structure for 4 base
%pair, HotMotifData.A or .C or .G or .T, which contains a Mx5 cell that
%stores the trinucloetide sequence + mutation position number, and the
%mutation of the .X nucleotide to A, C, G, or T.
%
%  HotMotifData = getHotMotifData(MotifData)
%
%  INPUT
%    MotifData: output structure from collectMotifData
%    
%  OUTPUT
%    HotMotifData: structure containing:
%      Header: {'TriNucCode' 'A' 'C' 'G' 'T'};
%      A: Mut freq for A -> X
%      C: Mut freq for C -> X
%      G: Mut freq for G -> X
%      T: Mut freq for T -> X

function HotMotifData = getHotMotifData(MotifData)
%Reg expression, Nth position that mutates
HotMotifs = {'[AT]A'     2;  %WA - ADAR
             '[AT][AG]C' 3;  %WRC - AID
             'G[CT][AT]' 1;  %GYW - AID
             'T[AT]'     1}; %TW - ADAR

%Mark which trinuc belong to which mutation         
TriNuc = MotifData.RowHeader;
MotifLabel = zeros(size(MotifData.Data, 1), 4);
for j = 1:length(TriNuc)
    for k = 1:size(HotMotifs, 1)
        HotLoc = regexp2(TriNuc{j}, HotMotifs{k, 1}) + HotMotifs{k, 2} - 1;
        MotifLabel(j, HotLoc) = k;
    end
end

%Extract the NNN# format for the trinucleotide + mut. position
LinTriNuc = repmat(MotifData.RowHeader, 1, 3);
for j = 1:length(LinTriNuc)
    for k = 1:3
        LinTriNuc{j, k} = [LinTriNuc{j, k} num2str(k)];
    end
end

%Format into a Mx5 cell with TriNucCode, Amut, Cmut, Gmut, Tmut
HotMotifData.Header = {'TriNucCode', 'A', 'C', 'G', 'T'};
for n = 1:4
    MotifIdx = find(MotifLabel == n);
    NucData = cell(length(MotifIdx), 4);
    for j = 1:length(MotifIdx)
        NucData(j, :) = num2cell(MotifData.Data{MotifIdx(j)});
    end
    HotMotifData.(int2nt(n)) = cat(2, LinTriNuc(MotifIdx), NucData);
end
