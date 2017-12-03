%getCleanGeneName will remove the - and * and / characters in the gene name
%with _ so that it can be used a structure field names.
%
%  GeneNames = getCleanGeneName(GeneNames)
%
%  INPUT
%    GeneNames: Mx1 cell of gene names from the databas of gene names
%      returned by DB = getGeneDatabase(Species). GeneNames would be
%      something like DB.Vmap(:, 2).
%   
%  OUTPUT
%    GeneNames: same size as input, but removeing non alpha-numeric
%      characters with '_'.

function GeneNames = getCleanGeneName(GeneNames)
for k = 1:length(GeneNames)
    GeneName = strrep(GeneNames{k}, ' ', '');
    ColLoc = find(GeneName == ':');
    if ~isempty(ColLoc)
        GeneName = GeneName(ColLoc(1)+1:end);
    end
    
    GeneName = regexp(GeneName, '\|', 'split');
    GeneName = GeneName{1};
    if(GeneName(5) == '0')
        GeneName(5) = [];
    end
    
    NonAlphaNum = regexp(GeneName, '\W');
    GeneName(NonAlphaNum) = '_';
    GeneNames{k} = GeneName;
end

