%getProfileStruct will return a structure from DB where the field names are
%those of the gene names with non-alpanumeric characters converted to '_',
%and the values of each field containing a double 1xN matrix where N is
%equal to the length of each reference gene. 
%
%  Prof = getProfileStruct(DB, Letter)
%
%  INPUT
%    DB: gene database    
%    Letter: gene type (V, D, J, Vk, Jk, Vl, Jl)
%
%  OUTPUT
%    Prof: structure containing cleaned gene names (getCleanGeneName) and
%      1xN double matrix to store SHM data.

function Prof = getProfileStruct(DB, Letter)
GeneNameLoc = contains(DB.MapHeader, 'GeneName', 'IgnoreCase', true);
RefSeqLoc = contains(DB.MapHeader, 'Seq', 'IgnoreCase', true);

ValidLetters = {'V', 'D', 'J', 'Vk', 'Jk', 'Vl', 'Jl'};
ValidLoc = find(ismember(lower(ValidLetters), lower(Letter)));
if isempty(ValidLoc)
    error('%s: not a valid Letter.', mfilename);
end
Map = DB.([ValidLetters{ValidLoc} 'map']);

Names = getCleanGeneName(Map(:, GeneNameLoc));
Len  = cellfun(@length, Map(:, RefSeqLoc));
Value  = cellfun(@(x) zeros(1,x), num2cell(Len), 'unif', false);
Temp = [Names Value]';
Prof = struct(Temp{:});
