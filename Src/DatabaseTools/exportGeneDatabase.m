%exportGeneDatabase will save the database sequences as a fasta file
%
%  exportGeneDatabase(Species, Format, varargin)
%
%  INPUT
%    Species: the species name for the database to export
%    Format ['imgt' or 'none']: the 'imgt' will add the gaps to the
%      sequences to be same as the IMGT gapp notation. Default is 'none'
%      where no gaps are added.
%    varargin: the  param-value pair for filtering databases (see
%      filterGeneDatabase.m).
% 
%  OUTPUT
%    Will ask user where to save the database to via a user interface.
%
function exportGeneDatabase(Species, Format, varargin)

DB = getGeneDatabase(Species);
DB = filterGeneDatabase(DB, varargin{:});

if nargin < 2
    Format = 'none';
end

MapNames = fieldnames(DB);
MapNames = MapNames(endsWith(MapNames, 'map'));

M = getMapHeaderVar(DB.MapHeader);

for j = 1:length(MapNames)
    DelLoc = cellfun(@isempty, DB.(MapNames{j})(:, M.Seq));
    DB.(MapNames{j})(DelLoc, :) = [];
end

[FileName, FilePath] = uiputfile('*.fa', 'Save Database File As');

for j = 1:length(MapNames)
    Xmap = DB.(MapNames{j});
    if isempty(Xmap); continue; end
    Xname = regexp(Xmap{1, M.Gene}, 'IG[HKL][VDJ]', 'match');
    Header = cell(size(Xmap, 1), 1);
    for k  = 1:size(Xmap, 1)
        Header{k} = sprintf('%s|', Xmap{k, [M.Gene, M.Funct, M.Strain]});
        Header{k}(end) = [];
    end
    if strcmpi(Format, 'imgt') %Add the gaps
        Sequence = cellfun(@(x, y) addIMGTgaps(x, y), Xmap(:, M.Seq), Xmap(:, M.GapInfo), 'un', false);
    else
        Sequence = Xmap(:, M.Seq);
    end
    OutFileName = fullfile(FilePath, [FileName(1:find(FileName == '.', 1, 'last')) Xname{1} '.fa']);
    fastawrite(OutFileName, Header, Sequence);
end