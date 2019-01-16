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

function G = getGrpIdx(varargin)
if iscell(varargin{1}) && isstruct(varargin{2}) || iscell(varargin{2})
    Map = getVDJmapper(varargin{2});
    VDJdata = varargin{1};
    GrpNum = cell2mat(VDJdata(:, Map.GrpNum));
    Template = cell2mat(VDJdata(:, Map.Template));
elseif isnumeric(varargin{1}) && isnumeric(varargin{2}) && numel(varargin{1}) == numel(varargin{2})
    GrpNum = varargin{1};
    Template = varargin{2};
else
    error('%s: Inputs are wrong. Must either be VDJdata and VDJheader or Map, or GrpNum and Template in Nx1 vectors.', mfilename);
end
%function G = getGrpIdx(VDJdata, VDJheader, SizeFilter)
%Map = getVDJmapper(VDJheader);

if nargin < 3 || isempty(varargin{3})
    SizeFilter = 'AC';
else
    SizeFilter = varargin{3};
end

%GrpNum = cell2mat(VDJdata(:, Map.GrpNum));
if isnumeric(SizeFilter)
    UnqGrpNum = SizeFilter; %It's not really SizeFilter, but the Group # you want to get
    Idx = cell(numel(UnqGrpNum), 1);
    for j = 1:numel(Idx)
        Idx{j} = find(GrpNum == UnqGrpNum(j));
    end
else
    [UnqGrpNum, ~, ~, Idx] = unique2(GrpNum);
    TransLoc = cellfun('size', Idx, 1) > 1; %transpot `1xM to Mx1 index, just in case
    Idx(TransLoc) = cellfun(@(x) x.', Idx(TransLoc), 'un', 0); %make sure it's 1xN vector to be able to do [G(1:N).Idx];
end

if isempty(UnqGrpNum)
    G = struct('GrpNum', [], 'Idx', [], 'Template', [], 'Size', []);
    return
end

G(1:length(UnqGrpNum)) = struct('GrpNum', [], 'Idx', [], 'Template', [], 'Size', []);
for y = 1:length(UnqGrpNum)
    G(y).GrpNum = UnqGrpNum(y);
    G(y).Idx = Idx{y}; 
    G(y).Template = sum(Template(G(y).Idx)); %sum(cell2mat(VDJdata(G(y).Idx, Map.Template)));
    G(y).Size = length(G(y).Idx);
end
if isnumeric(SizeFilter); return; end

SizeFilterStr = regexpi(SizeFilter, '[a-z]+', 'match', 'once');
N = convStr2NumMEX(SizeFilter);
switch upper(SizeFilterStr)
    case {'IND', 'CLONE'}
        G = cell2struct([num2cell(GrpNum), num2cell([1:numel(GrpNum)]'), num2cell(Template), num2cell(ones(numel(GrpNum), 1))], {'GrpNum', 'Idx', 'Template', 'Size'}, 2); %#ok<NBRAK>
        %G = cell2struct([VDJdata(:, Map.GrpNum) num2cell([1:size(VDJdata, 1)]') VDJdata(:, Map.Template) num2cell(ones(size(VDJdata, 1), 1))], {'GrpNum', 'Idx', 'Template', 'Size'}, 2); %#ok<NBRAK>
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