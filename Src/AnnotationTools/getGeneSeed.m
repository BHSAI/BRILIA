%getGeneSeed will look through the database V or J genes, and extract N
%number of nts to the left or right of the C or W codon set, including the
%104C and 118W codon set.
%
%  SeedSeq = getGeneSeed(DB, X, Nleft, Nright, Alphabet)
%
%  INPUT
%    DB: gene database structure(getCurrentDatabase.m)
%    X ['V' 'Vk' 'Vl' 'J' 'Jk' 'Jl']: Specificy what database gene is used.
%    Nleft: number of NT or AA left of anchor
%    Nright: number of NT or AA right of anchor
%    Alphabet ['nt' or 'aa']: Choose NT or AA letters
%
%  OUTPUT
%    SeedSeq: cell array of NT or AA sequences that can be used for
%      alignment purposes. 
%
%  NOTE
%    Genes without a solid anchor point C or W are REMOVED from SeedSeq
%    output! This is to prevent confusion pinpointing where the seed
%    sequence marks the position of the C or W.

function SeedSeq = getGeneSeed(DB, X, Nleft, Nright, Alphabet)
SeedSeq = {''};
if nargout < 5
    Alphabet = 'nt';
end

X = strsplit(strrep(X, ' ', ''), ',');
if any(startsWith(X, 'V')) > 1 && any(startsWith(X, 'J')) %Can't have both J and V
    error('%s: X input cannot have both V and J genes.', mfilename);
end

Fields = fieldnames(DB);
Fields(~endsWith(Fields, 'map')) = [];
ValidX = strrep(Fields, 'map', '');
GetIdx = find(endsWith(ValidX, X, 'ignorecase', true));
if isempty(GetIdx)
    error('%s: 2nd input X must be a comma-sperated list string containing [%s].', mfilename, makeStrPattern(strrep(Fields, 'map', ''), ','))
end

if length(GetIdx) > 1
    Xmap = cell(1, length(GetIdx));
    for j = 1:length(GetIdx)
        Xmap{j} = DB.(Fields{GetIdx(j)});
    end
    Xmap = vertcat(Xmap{:});
else
    Xmap = DB.(Fields{GetIdx});
end
if isempty(Xmap); return; end

%Remove all sequences without conserved residue or empty sequences
M = getMapHeaderVar(DB.MapHeader);
DelLoc = cellfun(@(x) x < 3, Xmap(:, M.Anchor)) | ...
         cellfun(@isempty, Xmap(:, M.Seq)) | ...
         ~strcmpi(Xmap(:, M.Funct), 'F');
Xmap(DelLoc,:) = [];
if isempty(Xmap); return; end

%Determine how many N number of nts to keep
if strcmpi(Alphabet, 'aa')
    Nleft = 3*Nleft;
    Nright = 3*Nright + 2; % The 2 comes from the anchor codon 2nd-3rd nt.
end

%Trim all sequence around the anchor point
TempSeedSeq = cell(size(Xmap, 1), 1);
if startsWith(X, 'V')
    TempSeedSeq = cellfun(@(x, y) padtrimSeq(x, length(x) - y + 1, Nleft, Nright), Xmap(:, M.Seq), Xmap(:, M.Anchor), 'un', false);
elseif startsWith(X, 'J')
    TempSeedSeq = cellfun(@(x, y) padtrimSeq(x, y, Nleft, Nright), Xmap(:, M.Seq), Xmap(:, M.Anchor), 'un', false);
end

if strcmpi(Alphabet, 'aa')
    TempSeedSeq = convNT2AA(TempSeedSeq, 'ACGTOnly', false);
end

TempSeedSeq = unique(TempSeedSeq);

%Determine similar by 1 sequences
[Ham, ~, ~, Penalty] = calcPairDistMEX(TempSeedSeq);
Ham(Penalty > 0) = Inf;
AncMap = calcAncMap(Ham);
ClustNum = findTreeClust(AncMap);

SeedSeq =  cell(max(ClustNum), 1);
for k = 1:max(ClustNum)
    ClustIdx = find(ClustNum == k);
    NTcount = countNT(TempSeedSeq(ClustIdx));
    MaxLoc = max(NTcount(1:4, :)) == length(ClustIdx);
    Seq = TempSeedSeq{ClustIdx(1)};
    Seq(~MaxLoc) = 'N';
    SeedSeq{k} = Seq;
end
