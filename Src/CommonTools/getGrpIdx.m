%getGrpIdx will find the indices of all the groups in the VDJdata set,
%returned as a structure.
%
%  G = getGrpIdx(VDJdata, VDJheader, SizeFilter)
%
%  G = getGrpIdx(VDJdata, VDJheader, GetGrpNum)
%
%  INPUT
%    VDJdata: BRILIA main data cell
%    VDJheader: BRILIA header cell
%    GetGrpNum: group number(s) to extract
%    SizeFilter: string code for filtering the based on clonotype sizes
%      'AC' - all clonotypes
%      'BC' - branched clonotypes with >= 2 unique sequences per clonotype
%      'BCN' - branched clonotypes with >= N unique sequences per clonotype
%      'TOPN' - top N clonotypes based on total template per clonotype
%      'BOTN' - top N clonotpyes based on total template per clonotype
%      'IND' - treat each sequence as its own individual group
%
%  OUTPUT
%    G: Nx1 structure storing the GrpNum and Idx
%      .GrpNum - Group number in VDJdata
%      .Idx - Index in VDJdata for the group
%      .Template - clonotype total template count
%      .Size - # of unique sequence per group (same as length(G(n).Idx))

function G = getGrpIdx(VDJdata, VDJheader, SizeFilter)
Map = getVDJmapper(VDJheader);

if nargin < 3 || isempty(SizeFilter)
    SizeFilter = 'AC';
end

GrpNum = cell2mat(VDJdata(:, Map.GrpNum));
if isnumeric(SizeFilter)
    UnqGrpNum = SizeFilter;
else
    UnqGrpNum = unique(GrpNum);
end

if isempty(UnqGrpNum)
    G = struct('GrpNum', [], 'Idx', [], 'Template', [], 'Size', []);
    return
end

G(1:length(UnqGrpNum)) = struct('GrpNum', [], 'Idx', [], 'Template', [], 'Size', []);
for y = 1:length(UnqGrpNum)
    G(y).GrpNum = UnqGrpNum(y);
    G(y).Idx = find(UnqGrpNum(y) == GrpNum)'; %make sure it's 1xN vector
    G(y).Template = sum(cell2mat(VDJdata(G(y).Idx, Map.Template)));
    G(y).Size = length(G(y).Idx);
end
if isnumeric(SizeFilter); return; end

SizeFilterStr = regexpi(SizeFilter, '[a-z]+', 'match', 'once');
N = convStr2NumMEX(SizeFilter);
switch upper(SizeFilterStr)
    case 'IND'
        G = cell2struct([VDJdata(:, Map.GrpNum) num2cell([1:size(VDJdata, 1)]') VDJdata(:, Map.Template) num2cell(ones(size(VDJdata, 1), 1))], {'GrpNum', 'Idx', 'Template', 'Size'}, 2); %#ok<NBRAK>
        return
    case 'AC'
        return
    case 'BC'
        if isempty(N); N = 2; end
        BCLoc = cellfun('length', {G.Idx}) >= N;
        G = G(BCLoc);
        return
    case 'SC'
        SCLoc = cellfun('length', {G.Idx}) == 1;
        G = G(SCLoc);
        return
    case 'TOP'
        assert(~isempty(N), '%s: TopN must have a number >= 0 for N. Ex: ''Top100''', mfilename);
        [~, GetIdx] = findTopN([G.Template], N);
        G = G(GetIdx);
        return
    case 'BOT'
        assert(~isempty(N), '%s: BotN must have a number >= 0 for N. Ex: ''Bot100''', mfilename);
        [~, GetIdx] = findBotN([G.Template], N);
        G = G(GetIdx);
        return
    otherwise
        error('%s: Unrecognized SizeFilter option "%s".', mfilename, SizeFilter);
end