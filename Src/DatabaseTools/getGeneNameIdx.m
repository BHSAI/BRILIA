%getGeneNameIdx will intake a cell of gene names and return a matrix of
%indices pointing to the nth item in a gene database structure for a
%particular gene set. This is used mainly to function without BRILIA's
%"xGeneNum" data column, which is used only by BRILIA and serves no other
%purpose. 
%
%  Idx = getGeneNameIdx(GeneName, DB, X)
%
%  INPUT
%    GeneName: char or Mx1 cell of char of gene names.
%    DB: the database structure from getGeneDatabase.m
%    X['V' 'D' 'J' 'Vk' 'Vl' 'Jk' 'Jl']: the type of genes. Ex: Vk = light
%      chain kappa V genes
%
%  NOTE
%    If GeneName has multiple options like "IGHV1-11*01|IGHV1-11*02", then
%    this will only return the index of the first gene name.
%
function Idx = getGeneNameIdx(GeneName, DB, X)
if ischar(GeneName)
    GeneName = {GeneName};
end
assert(min(size(GeneName)) == 1, '%s: GeneName input must be a 1xN or Mx1 cell.', mfilename)

F = fieldnames(DB);
F = F(endsWith(F, 'map'));
XType = strrep(F, 'map', '');
MapIdx = find(endsWith(XType, X, 'ignorecase', true), 1);
assert(~isempty(MapIdx), '%s: 3rd input X is not valid. Must be [%s].', mfilename, makeStrPattern(XType, ','))

[UnqGenes, ~, UnqIdx] = unique(regexprep(GeneName, '\|.*', ''));
M = getMapHeaderVar(DB.MapHeader);
[~, UnqGeneIdx, GeneIdx] = intersect(UnqGenes, DB.(F{MapIdx})(:, M.Gene));

Idx = zeros(size(GeneName));
for j = 1:length(UnqGeneIdx)
    Idx(UnqIdx == UnqGeneIdx(j)) =  GeneIdx(j);
end