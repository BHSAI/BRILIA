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
    StrainName = sprintf(StrPat, Xmap{:, M.StrainLoc});   
    StrainCell = unique(regexp(StrainName, ';', 'split'));
    AllStrainCell = [AllStrainCell(:); StrainCell(:)];
end
UnqStrain = unique(AllStrainCell);
if isempty(UnqStrain{1})
    UnqStrain(1) = [];
end

%Regroup strain names based on the first X condense letter.
if N > 0 && length(UnqStrain) > 1
    %Prepare the UnqStrain for doing pairwise hamming distance matching
    UnqStrainChar = char(UnqStrain);
    UnqStrainChar = UnqStrainChar(:, 1:N);
    UnqStrainChar(UnqStrainChar == 'X') = '@'; %Ham dist calc assumes x is wildcard, so go and change these X out
    UnqStrainChar(UnqStrainChar == 'x') = '@'; 
    UnqStrainChar(UnqStrainChar == ' ') = 'X'; %Assume empty spaces are wildcards
    
    %Perform pairwise hamming distance calc to determine word clustering
    PairDist = calcPairDist(UnqStrainChar, 'ham');
    PairDistMap = calcAncMap(PairDist);
    PairDistMap(PairDistMap(:, 3) > 0, 2) = 0; 
    ClustMap = findTreeClust(PairDistMap);
    
    %Cluster strains and subsets of strain into the same cell
    UnqStrainCell = cell(max(ClustMap(:, 2)), 1);
    for j = 1:length(UnqStrainCell)
        StrainIdx = find(ClustMap(:, 2) == j);
        RepPat = repmat('%s;', 1, length(StrainIdx));
        RepPat(end) = [];
        UnqStrainCell{j} = sprintf(RepPat, UnqStrain{StrainIdx});
    end
else
    UnqStrainCell = UnqStrain;
end
