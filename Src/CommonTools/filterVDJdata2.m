function VDJdata = filterVDJdata2(VDJdata, VDJheader, varargin)
P = inputParser;
addParameter(P, 'Logic', 'and', @(x) ismember(lower(x), {'and', 'or'}));
addParameter(P, 'GetGrpNum', 0, @isnumeric);
addParameter(P, 'GetSeqNum', 0, @isnumeric);
addParameter(P, 'hCDR3', 0, @ischar);
addParameter(P, 'lCDR3', 0, @ischar);
addParameter(P, 'GrpMinSize', 0, @isnumeric);
addParameter(P, 'GrpMaxSize', Inf, @isnumeric);
addParameter(P, 'hCDR3MinSize', 0, @isnumeric);
addParameter(P, 'hCDR3MaxSize', Inf, @isnumeric);
addParameter(P, 'lCDR3MinSize', 0, @isnumeric);
addParameter(P, 'lCDR3MaxSize', Inf, @isnumeric);
parse(P, varargin{:});
P = P.Results;
Logic = P.Logic;
GetSeqNum = P.GetSeqNum;
GetGrpNum = P.GetGrpNum;
hCDR3 = P.hCDR3;
lCDR3 = P.lCDR3;
GrpMinSize = P.GrpMinSize;
GrpMaxSize = P.GrpMaxSize;
hCDR3MinSize = P.hCDR3MinSize;
hCDR3MaxSize = P.hCDR3MaxSize;
lCDR3MinSize = P.lCDR3MinSize;
lCDR3MaxSize = P.lCDR3MaxSize;

Map = getVDJmapper(VDJheader);

j = 1;

GrpNum = cell2mat(VDJdata(:, Map.GrpNum));
UnqGrpNum = unique(GrpNum);
if GrpMinSize > 0 || GrpMaxSize < Inf
    TmpLocs = zeros(size(VDJdata, 1), 1, 'logical');
    for y = 1:length(UnqGrpNum)
        Idx = UnqGrpNum(y) == GrpNum;
        Size = sum(Idx);
        if Size >= GrpMinSize && Size <= GrpMaxSize
            TmpLocs(Idx) = 1;
        end
    end
    Locs{j} = TmpLocs;
    j = j+1;
end

if GetGrpNum ~= 0
    Locs{j} = cellfun(@(x) isequal(x, GetGrpNum), VDJdata(:, Map.GrpNum));
    j = j+1;
end

if GetSeqNum ~= 0
    Locs{j} = cellfun(@(x) isequal(x, GetSeqNum), VDJdata(:, Map.SeqNum));
    j = j+1;
end

if contains(Map.Chain, 'h', 'ignorecase', true)
    if hCDR3MinSize > 0 || hCDR3MaxSize < Inf
        Locs{j} = cellfun(@(x) x >= hCDR3MinSize && x <= hCDR3MaxSize, VDJdata(:, Map.hCDR3(2)));
        j = j+1;
    end
    if ~isempty(hCDR3)
        Locs{j} = cellfun(@(x) strcmpi(x, hCDR3), VDJdata(:, Map.hCDR3(1)));
        j = j+1;
    end
end

if contains(Map.Chain, 'l', 'ignorecase', true)
    if lCDR3MinSize > 0 || lCDR3MaxSize < Inf
        Locs{j} = cellfun(@(x) x >= lCDR3MinSize && x <= lCDR3MaxSize, VDJdata(:, Map.lCDR3(2)));
        j = j+1;
    end
    if ~isempty(lCDR3)
        Locs{j} = cellfun(@(x) strcmpi(x, lCDR3), VDJdata(:, Map.lCDR3(1)));
        j = j+1;
    end
end

if j == 1
    return
elseif j == 2
    GetLocs = Locs{1};
else
    if strcmpi(Logic, 'and')
        GetLocs = and(Locs{:});
    else
        GetLocs = or(Locs{:});
    end
end

GetUnqGrp = unique(GrpNum(GetLocs));
for y = 1:length(GetUnqGrp)
    GetLocs = GetLocs | (GrpNum == GetUnqGrp(y));
end

VDJdata = VDJdata(GetLocs, :);