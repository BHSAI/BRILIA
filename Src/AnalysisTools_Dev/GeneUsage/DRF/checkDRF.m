%checkDRF will check the database D genes of a species to identify which
%reading frame are valid and what sequences each reading frame translates.
%
%  Summary = checkDRF(Species)
%
%  INPUT
%    Species: species of the database to upload
%
%  OUTPUT
%    Summary: structure of results
%      .AA - Mx4 cell, where col 2:4 are the RF1, 2, 3 amino acid
%        translation, and col 1 is the D gene name.
%      .NoStop - Mx4 cell, where col 1 = D gene name,  col 2:4 are 1 or 0
%        if there is no stop (1) or stop (0) codons.
%      .Dmap - MxN cell matrix of the D genes used
function Summary = checkDRF(Species)

DB = getGeneDatabase(Species);
Dmap = DB.Dmap;
[~, SortIdx] = sort2(Dmap(:, 2));
Dmap = Dmap(SortIdx, :);

AA = cell(size(Dmap, 1), 4);
AA(:, 1) = Dmap(:, 2);
for k = 1:size(Dmap, 1)
    AA (k, 2:end) = nt2aa(Dmap{k, 1}, 'acgtonly', false, 'frame', 'all');
end

NoStop = [AA(:, 1) num2cell(cellfun(@(x) ~contains(x, '*'), AA(:, 2:end)))];
Summary = struct;
Summary.AA = AA;
Summary.NoStop = NoStop;
Summary.Dmap = Dmap;
