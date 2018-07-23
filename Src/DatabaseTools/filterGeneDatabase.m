%filterRefGene will process the database DB to retain only those of
%interest, as defined by the filter criteria.
%
%  DB = filterRefGene(DB, Param, Value)
%  
%  INPUT
%    DB: database structure from getCurrentDatabase.m
%
%     Param       Value (defaults = *)     Details
%     ----------- ------------------------ --------------------------------
%     Strain      * all                    For mouse only. Use all strains. 
%                   c57bl                  C57BL strains, include C57BL6, C57BL/J
%                   balb                   BALB strains, such as BALB/C
%                   ExactName              NOTE: Some databases are incomplete for certain strains
%     Dgene       * all                    Foward and inverse are okay
%                   fwd                    Foward only
%                   inv                    Inverse only
%                   rev                    Inverse only (same as 'inv'. Will be phased out)
%     Vgene       * f                      Functional only
%                   p                      Psueudo genes only
%                   orf                    Open reading frame only
%                   all                    All of the above
%
%  OUTPUT
%    DB: the filtered version of the input DB, where unused gene sequences
%      are set to empty ''. This prevents the shifting of Map numbering
%      schemes required by the software to perform.
 
function [DB, FiltOption] = filterGeneDatabase(DB, varargin)
P = inputParser;
addParameter(P, 'Strain', 'all', @ischar);
addParameter(P, 'Dgene',  'all', @(x) ischar(x) && ismember(lower(x), {'all', 'fwd', 'inv'}));
addParameter(P, 'Vgene',  'f',   @(x) ischar(x) && all(ismember(strsplit(lower(x), ','), {'all', 'f', 'p', 'orf'})));
parse(P, varargin{:});
P = P.Results;

Fields = fieldnames(DB);
Fields = Fields(endsWith(Fields, 'map'));
M = getMapHeaderVar(DB.MapHeader);

FiltOption = struct('Vgene', 'f', 'Dgene', 'all', 'Strain', 'all');

%If dealing with mouse database, select by mouse strain.
if contains(DB.FilePath, fullfile('Databases', 'Mouse'), 'ignorecase', true)
    UnqStrain = [{'all'}; getUnqStrain(DB, 4)];
    if isempty(P.Strain) %Ask users to select desired strain
        StrainNum = chooseFromList(UnqStrain, 'Attempt', 5, 'Default', 1, 'Message', 'What mouse strain is it?');
    else %Determine which unq strain matches with user input
        StrainNum = find(contains(UnqStrain, P.Strain, 'ignorecase', true));
        if isempty(StrainNum)
            fprintf('Could not find a match for "%s".\n', P.Strain);
            StrainNum = chooseFromList(UnqStrain, 'Attempt', 5, 'Default', 1, 'Message', 'What mouse strain is it?');
        elseif length(StrainNum) > 1 %Use must decide which one to pick
            fprintf('Found multiple matches for "%s".\n', P.Strain);
            StrainNumT = chooseFromList(UnqStrain(StrainNum), 'Attempt', 5, 'Default', 1, 'Message', 'What mouse strain is it?');
            StrainNum = StrainNum(StrainNumT);
        end
    end
    
    %Filter database by setting the sequence to empty for non-matching strains
    if StrainNum > 1 
        SearchFor = strrep(UnqStrain{StrainNum}, ';', '|');
        for j = 1:length(Fields)
            DelLoc = cellfun(@isempty, regexpi(DB.(Fields{j})(:, M.Strain), SearchFor, 'once'));
            DB.(Fields{j})(DelLoc, M.Seq) = {''};
        end
    end
    FiltOption.Strain = UnqStrain{StrainNum};  %Save the final filter value
end

%==========================================================================
%Determine if you want to search Dinverse too
DgeneList = {'all'; 'fwd'; 'inv'};
if isempty(P.Dgene) %Ask user to determine which direction D to search
    DgeneNum = chooseFromList(DgeneList, 'Attempt', 5, 'Default', 1, 'Message', 'Which D gene directions?');
else %Search of user selection
    DgeneNum = find(ismember(DgeneList, P.Dgene));
    if isempty(DgeneNum)
        DgeneNum = chooseFromList(DgeneList, 'Attempt', 5, 'Default', 1, 'Message', 'Which D gene directions?');
    end
end
switch DgeneNum
    case 2 %Keep only Dfwd by deleting inverse
        DB.Dmap( startsWith(DB.Dmap(:, M.Gene), 'r'), M.Seq) = {''};
    case 3 %Keep only Dinv by deleting forward
        DB.Dmap(~startsWith(DB.Dmap(:, M.Gene), 'r'), M.Seq) = {''};
end
FiltOption.Dgene = DgeneList{DgeneNum};

%==========================================================================
%Determine if you want to include pseudo or orf sequences
VgeneList = {'all'; 'f'; 'p'; 'orf'};
if isempty(P.Vgene) %Ask user to determine which direction D to search
    VgeneNum = chooseFromList(VgeneList, 'Attempt', 5, 'Default', 2, 'MultiSelect', 'on', 'Message', 'Which V gene functionality? Multiple allowed (ex: 2,3).');
else %Search of user selection
    VgeneNum = find(ismember(VgeneList, P.Vgene));
    if isempty(VgeneNum)
        VgeneNum = chooseFromList(VgeneList, 'Attempt', 5, 'Default', 2, 'MultiSelect', 'on', 'Message', 'Which V gene functionality? Multiple allowed (ex: 2,3).');
    end
end
if min(VgeneNum) > 1 %You have to filter
    Idx = find(startsWith(Fields, 'V', 'ignorecase', true));
    for j = 1:length(Idx)
        Funct = regexpi(DB.(Fields{Idx(j)})(:, M.Funct), '\w+', 'match');
        Funct = cellfun(@(x) lower(x{1}), Funct, 'un', 0);
        DelLoc = ~ismember(Funct, VgeneList(VgeneNum));
        DB.(Fields{Idx(j)})(DelLoc, M.Seq) = {''};
    end
end
FiltOption.Vgene = makeStrPattern(VgeneList(VgeneNum), ',');