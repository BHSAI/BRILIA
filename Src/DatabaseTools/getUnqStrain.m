%getUnqStrain will look at the VDJ database to look for unique mouse
%strains. Returns a cell with  unique strains.
%
%  UnqStrainCell = getUnqStrain(DB)
%
%  UnqStrainCell = getUnqStrain(DB, N)
%
%  INPUT
%    DB: database structure from getGeneDatabase
%    N: the first N number of letters of the strain name to use for
%      condensing the strings. If N = 0, will only return unique, full name
%      strains. If N = 4, will match strings up to N letter and group them
%      together. 
%
%  OUTPUT
%    UnqStrianCell: a Mx1 cell of string of unique strains, group based on
%      how many N letters of the strain name are the same.
%
%  EXAMPLE
%    DB = getGeneDatabase('Mouse')
%    UnqStrainCell = getUnqStrain(DB, 0)
%    UnqStrainCell = 
%         '129/Sv'
%         'A/J'
%         'AKR'
%         'BALB.K'
%         'BALB/b'
%         'BALB/c'
%         'C3H'
%         'C57BL'
%         'C57BL/10'
%         'C57BL/6'
%         'C57BL/6J'
%         ...
%
%    UnqStrainCell = getUnqStrain(DB, 4)
%    UnqStrainCell = 
%         '129/Sv'
%         'A/J'
%         'AKR'
%         'BALB.K;BALB/b;BALB/c'
%         'C3H'
%         'C57BL;C57BL/10;C57BL/6;C57BL/6J'
%         ...

function UnqStrainCell = getUnqStrain(DB, N)
%Determine map fields only
M = getMapHeaderVar(DB.MapHeader);
Fields = fieldnames(DB);
MapLoc = findCell(Fields, 'map', 'MatchWord', 'partial');
Fields = Fields(MapLoc);

%Collect all unique strain names
AllStrainCell = {};
for j = 1:length(Fields)
    Xmap = DB.(Fields{j});
    if isempty(Xmap)
        continue; 
    end
    StrPat = repmat('%s;', 1, size(Xmap, 1));
    StrPat(end) = [];
    StrainName = sprintf(StrPat, Xmap{:, M.Strain});   
    StrainCell = unique(regexp(StrainName, ';', 'split'));
    AllStrainCell = [AllStrainCell(:); StrainCell(:)];
end
UnqStrain = unique(AllStrainCell);
if isempty(UnqStrain{1})
    UnqStrain(1) = [];
end

%Regroup strain names based on the first X condense letter.
if N > 0 && length(UnqStrain) > 1
    UnqStrainN = cell(length(UnqStrain), 1);
    for j = 1:length(UnqStrain)
        if length(UnqStrain{j}) >  N
            UnqStrainPre = UnqStrain{j}(1:N);
        else
            UnqStrainPre = UnqStrain{j};
        end
        UnqStrainPre = regexp(UnqStrainPre, '\w+', 'match');
        UnqStrainN{j} = UnqStrainPre{1};
    end
    
    [~, ~, UnqIdx] = unique(UnqStrainN);
    UnqStrainCell = cell(max(UnqIdx), 1);
    for j = 1:length(UnqStrainCell)
        Idx = UnqIdx == j;
        TempUnqStr = sprintf('%s;', UnqStrain{Idx});
        UnqStrainCell{j} = TempUnqStr(1:end-1);
    end
else
    UnqStrainCell = UnqStrain;
end
