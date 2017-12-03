%getGeneUsage will plot a VDJ gene usage in either a scatter plot or
%histogram type setting. The gene can either be done at the family number
%or full gene name.
%
%  INPUT
%    FileNames: BRILIA output file name to process
%    DetailLevel ['all' or 'family']: all gene name or by gene family #
%
function GeneUsage = getGeneUsage(FileNames, DetailLevel)
if nargin == 0
    FileNames = '';
end

if nargin < 2
    DetailLevel = 'all';
end

FileNames = getBriliaFiles(FileNames);
if ischar(DetailLevel)
    DetailType = {'all', 'family'};
    DetailLoc = find(startsWith(DetailType, DetailLevel, 'IgnoreCase', true));
    assert(~isempty(DetailLoc), '%s: Could not determine DetailLevel "%s". Must be "%s" or "%s".', mfilename, DetailLevel, DetailType{:});
    DetailLevel = DetailType(DetailLoc(1));
end

for f = 1:length(FileNames)
    [VDJdata, VDJheader] = openSeqData(FileNames{f});
    Map = getVDJmapper(VDJheader);
    switch Map.Chain
        case 'H'
            GeneLoc = Map.hGeneName;
        case 'L'
            GeneLoc = Map.lGeneName;
        case 'HL'
            GeneLoc = [Map.hGeneName Map.lGeneName];
    end
    NumGenes = length(GeneLoc);
    
    GrpNum = cell2mat(VDJdata(:, Map.GrpNum));
    [~, UnqIdx, ~] = unique(GrpNum);
    
    if f == 1
        UnqVDJ = VDJdata(UnqIdx, GeneLoc);
    else
        UnqVDJ = cat(1, UnqVDJ, VDJdata(UnqIdx, GeneLoc));
    end
end
NumClones = size(UnqVDJ, 1);

%Count up the occurence of each V/D/J genes, counting ambiguious hits as
%1/N frequency, where N is the "number of VDJ+VJ permutations allowed in
%the overall anotation of that clone. (Example: For V1, D1|D2|D3, J1|J2, N
%would be 6.
UnqGeneNames = cell(1, NumGenes);
for n = 1:NumGenes
    AllNames = sprintf('%s,', UnqVDJ{:, n});
    UnqGeneNames{n} = unique(regexp(AllNames(1:end-1), ',|\|', 'split'));
    if strcmpi(DetailLevel, 'family')
        Family = unique(parseGeneName(UnqGeneNames{n}));
        FamilyNum = regexp(Family, '\d+', 'match');
        [~, SortIdx] = sort(str2double([FamilyNum{:}]));
        UnqGeneNames{n} = Family(SortIdx);
    end
end

GeneCombCts = zeros(cellfun(@length, UnqGeneNames));
GeneIdx = cell(NumGenes, 1);
for j = 1:NumClones
    for n = 1:NumGenes
        GeneName = regexp(UnqVDJ{j, n}, '\|', 'split');
        if strcmpi(DetailLevel, 'family')
            GeneName = unique(parseGeneName(GeneName));
        end
        GeneIdx{n} = find(ismember(UnqGeneNames{n}, GeneName));
    end
    GeneCombVec = combvec2(GeneIdx{:});
    GeneCombCell = cell(NumGenes, 1);
    for k = 1:NumGenes
        GeneCombCell{k} = GeneCombVec(k, :);
    end
    Idx = sub2ind(size(GeneCombCts), GeneCombCell{:});
    GeneCombCts(Idx) = 1/length(Idx);
end

GeneUsage.Names = UnqGeneNames;
GeneUsage.Count = GeneCombCts;