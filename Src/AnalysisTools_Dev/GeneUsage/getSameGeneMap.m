%getSameGeneMap will return a Mx2 cell of gene names, where the 1st
%column is the gene name and the 2nd column is the representative gene name
%to use due to 100% identiical nt sequence.
%
%  [MapAll2Unq, MapUnq2All] = getSameGeneMap(Species, Gene)
%
%  INPUT
%    Species: species of the database to use
%    Gene ['V' 'D' 'J' 'Vk' 'Jk' 'Vl' 'Jl']: gene to process
%
%  OUTPUT
%    MapAll2Unq: Mx2 cell where col1 = every gene and col2 = unique
%       repsentative Gene (U is less than M)
%    MapUnq2All: Ux2 cell where col1 = representative unqiue gene and col2
%       = cell of possible unique gene (U is less than M)

function [MapAll2Unq, MapUnq2All] = getSameGeneMap(Species, Gene)
DB = getGeneDatabase(Species);

ValidGene = {'V' 'D' 'J' 'Vk' 'Jk' 'Vl' 'Jl'};
[~, ~, ValidIdx] = intersect(lower(Gene), lower(ValidGene));
if isempty(ValidIdx)
    error('%s: Invalid Gene given "%s". Should be any of: %s.', mfilename, Gene, makeStrPattern(ValidGene, ',')); 
end
Gene = ValidGene{ValidIdx};
GeneMap = DB.([Gene 'map']);

NameIdx = findCell(DB.MapHeader, 'GeneName');
SeqIdx = findCell(DB.MapHeader, 'Seq');
[~, UnqNameIdx, UnqIdx] = unique(GeneMap(:, SeqIdx), 'stable');
UnqName = GeneMap(UnqNameIdx, NameIdx);

MapAll2Unq = cell(size(GeneMap, 1), 2);
MapAll2Unq(:, 1) = GeneMap(:, NameIdx);
MapAll2Unq(:, 2) = UnqName(UnqIdx);% 

if nargout == 2
    MapUnq2All = cell(max(UnqIdx), 2);
    MapUnq2All(:, 1) = GeneMap(UnqNameIdx, NameIdx);
    for k = 1:max(UnqIdx)
        MapUnq2All{k, 2} = GeneMap(k == UnqIdx, NameIdx);
    end
end