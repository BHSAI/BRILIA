%countGeneFreq will count the frequence of the VDJ/VJ genes in the data set.
%
%  INPUT
%    VDJdata: BRILIA data cell
%    VDJdata: BRILIA data header cell
%
%  OUTPUT
%    GeneFreq: structure storing the frequency results
%      .Fam.[hV/hD/hJ/lV/lJ] - Nx2 cell of unique gene family frequency
%      .All.[hV/hD/hJ/lV/lJ] - Nx2 cell of unique gene name frequency

function GeneFreq = countGeneFreq(VDJdata, VDJheader)
Map = getVDJmapper(VDJheader);
G = getGrpIdx(VDJdata, VDJheader);
GeneFreq = struct('Fam', {}, 'All', {});

for c = 1:length(Map.Chain)
    FieldName = [lower(Map.Chain(c)) 'GeneName'];
    
    for n = 1:length(Map.(FieldName))
        GeneData = cell(length(G), 2); %GeneFam, GeneName (1st one only)
        for y = 1:length(G)
            GeneName = regexp(VDJdata{G(y).Idx(1), Map.(FieldName)(n)}, '\|', 'split');
            GeneFam = parseGeneName(GeneName{1});
            GeneData(y, :) = {GeneFam{1} GeneName{1}};
        end
        
        [~, SortIdx] = sortByFamNum(GeneData(:, 1));
        GeneData = GeneData(SortIdx, :);

        [UnqGeneFam,  ~, UnqFamIdx] = unique(GeneData(:, 1), 'stable');
        FamFreq = zeros(length(UnqGeneFam), 1);
        for k = 1:length(UnqGeneFam)
            FamFreq(k) = sum(k == UnqFamIdx);
        end

        [UnqGeneName, ~, UnqNameIdx] = unique(GeneData(:, 2), 'stable');
        NameFreq = zeros(length(UnqGeneName), 1);
        for k = 1:length(UnqGeneName)
            NameFreq(k) = sum(k == UnqNameIdx);
        end
        
        if strcmpi(UnqGeneName{1}(1), 'r') %in case reverse rIGHD, need 5th for D
            M = 5;
        else
            M = 4;
        end
        GeneFreq(1).Fam.([lower(Map.Chain(c)) UnqGeneName{1}(M)]) = [UnqGeneFam  num2cell(FamFreq)];
        GeneFreq(1).All.([lower(Map.Chain(c)) UnqGeneName{1}(M)]) = [UnqGeneName num2cell(NameFreq)];
    end
end