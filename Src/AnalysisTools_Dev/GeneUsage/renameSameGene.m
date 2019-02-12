%renameSameGene will rename genes that have the same nt sequence into a
%unique representative gene.
%
%  GeneName = renameSameGene(GeneName, Species)
%
%  S = renameSameGene(S, Species)
%
%  INPUT
%    S: structure of VDJdata and VDJheader
%    GeneName: char or cell of string of gene names
%    Species: species of the database to use
%
%  OUPUT
%    GeneName: char or cell of string of gene names where degenerate gene
%      names are relabeled
%    S: same as input S, but all VDJ heavy and light chain names are
%      altered to unique representative names

function Out = renameSameGene(varargin)
if nargin < 1 || nargin >= 3
    error('%s: Need 1 or 2 inputs.', mfilename);
end
if isstruct(varargin{1}) 
    Out = varargin{1};
    if nargin == 2 && ~isfield(varargin{1}, 'Species') && ischar(varargin{2})
        for f = 1:length(Out)
            Out(f).Species = varargin{2};
        end
    end
    if ~all(isfield(Out, {'VDJdata', 'VDJheader', 'Species'}))
        error('%s: When using structure input, there must be VDJdata, VDJheader, and Species fields.', mfilename);
    end
    for f = 1:length(Out)
        Map = getVDJmapper(Out(f).VDJheader);
        GeneIdx = vertcat(Map.hGeneName(:), Map.lGeneName(:));
        GeneIdx(GeneIdx == 0) = [];
        for k = 1:length(GeneIdx)
            Out(f).VDJdata(:, GeneIdx(k)) = renameSameGene_T(Out(f).VDJdata(:, GeneIdx(k)), Out(f).Species);
        end
    end
elseif nargin == 2 && iscell(varargin{1}) && ischar(varargin{2}) 
    Out = renameSameGene_T(varargin{1}, varargin{2});
else
    error('%s: Recheck the inputs.', mfilename);
end

function GeneName = renameSameGene_T(GeneName, Species)
GeneType = {'IGHV' 'IGHD' 'IGHJ' 'IGKV' 'IGKJ' 'IGLV' 'IGLJ';
            'V'    'D'    'J'    'Vk'   'Jk'   'Vl'   'Jl' };
Gene = '';
for k = 1:size(GeneType, 2)
    if contains(GeneName, GeneType{1, k})
        Gene = GeneType{2, k};
        GeneName = renameGeneName_T2(GeneName, Species, Gene);
    end
end
if isempty(Gene)
    error('%s: Could not determine the gene type from the GeneName given.', mfilename);
end

function GeneName = renameGeneName_T2(GeneName, Species, Gene)
IsChar = ischar(GeneName);
if IsChar
    GeneName = {GeneName};
end

MapAll2Unq = getSameGeneMap(Species, Gene);
if Gene == 'D' %fwd and rev same map are same. so keep one.
    RevLoc = cellfun(@(x) x(1) == 'r', MapAll2Unq(:, 1));
    MapAll2Unq(RevLoc, :) = [];
end

for k = 1:size(MapAll2Unq, 1)
    GeneName = strrep(GeneName, MapAll2Unq{k, 1}, MapAll2Unq{k, 2});
end

MultLoc = contains(GeneName, '|');
GeneName(MultLoc) = cellfun(@(x) unique(regexp(x, '\||\,', 'split')), GeneName(MultLoc), 'unif', false);
GeneName(MultLoc) = cellfun(@(x) x{1}, GeneName(MultLoc), 'unif', false);

if IsChar
    GeneName = GeneName{1};
end