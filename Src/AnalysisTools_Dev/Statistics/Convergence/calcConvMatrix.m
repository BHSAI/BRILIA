%calcConvMatrix will collect a 1xN cell of Mx2 cell of [Entity, Weight],
%where Entity is a sequence (case-senstive) and Weight is the weight of
%that sequence.
%
%  [Heatmap, Lookup, Template] = calcConvMatrix(Data1, ..., Data3)
%
%  INPUT
%    DataN: Mx2 cell where col 1 is the sequence entity, col 2 is the count
%    of that sequence
%
%  OUTPUT
%    Heatmap: structure of various NxN matrix of repertoire comparison
%    metrics/premetrics
%      .Ovlerlap      # entities overlap / total # entities
%      .WtOverlap     # weighted entities overlap / total # weighted entities
%      .BhattDist     Bhattacharya distance of weights of overlapping entities
%      .KLDivergence  Kullback-Leibner diverence of weights of overlapping entities
%
%  EXAMPLE
%    Mat{1} = {'a' 1; 'b' 2; 'c' 1; 'd' 3};
%    Mat{2} = {'a' 1; 'b' 2; 'e' 1; 'f' 7};
%    Mat{3} = {'c' 1; 'd' 3; 'e' 1; 'f' 7}; 
%    Conv = calcConvMatrix(Mat{:});
%    Conv.SeqOvlerlap = 
%     0         0.5000    0.5000
%     0.5000         0    0.5000
%     0.5000    0.5000         0
%
%    Conv.CellOverlap = 
%          0    0.4286    0.5714
%     0.2727         0    0.7273
%     0.3333    0.6667         0
%
%    Conv.BhattDist = 
%      0     0     0
%      0     0     0
%      0     0     0
%
%    Conv.KLDivergence = 
%      0     0     0
%      0     0     0
%      0     0     0
function Out = calcConvMatrix(varargin)
Data = cell(1, nargin);
Weight = cell( 1, nargin);
for j = 1:nargin
    Data{j} = varargin{j}(:, 1);
    if size(varargin{j}, 2) == 2
        Weight{j} = cell2mat(varargin{j}(:, 2));
    else
        Weight{j} = ones(size(Data{j}));
    end
end
[UnqSeq, UnqIdx] = getUnqSeq(Data);
Lookup = getLookup(Data, UnqSeq, UnqIdx);
LookupWeight = lookupData(Data, Lookup, Weight);
Out = calcConv(Data, Weight, Lookup, LookupWeight);

function [UnqSeq, UnqIdx] = getUnqSeq(Data)
N = numel(Data);
UnqSeq = cell(N, 1);
UnqIdx = cell(N, 1);
for f = 1:N
    [UnqSeq{f}, ~, ~, UnqIdx{f}] = unique2(Data{f});
end

%Create a NxN cell of indices of overlapping entities
function Lookup = getLookup(Data, UnqSeq, UnqIdx)
N = numel(Data);
Lookup = cell(N);
fprintf('%s: Finding overlaps: %d...', mfilename, 1);
for r = 2:N
    fprintf('%d...', r);
    for c = 1:r-1
        [~, IA, IB] = intersect(UnqSeq{r}, UnqSeq{c});
        Lookup{r, c} = UnqIdx{r}(IA);
        Lookup{c, r} = UnqIdx{c}(IB);
    end
end
fprintf('\n');

%Generates a NxN cell array, replacing index in Lookup with the value in T
function Out = lookupData(Data, Lookup, T)
N = numel(Data);
Out = cell(N);
for r = 1:N
    Temp = T{r};
    for c = 1:N
        if r == c; continue; end
        Out{r, c} = cellfun(@(x) sum(Temp(x)), Lookup{r, c});
    end
end
    
function Out = calcConv(Data, Weight, Lookup, LookupWeight)
N = numel(Data);
Out.SeqOverlap = zeros(N);
Out.CellOverlap = zeros(N);
Out.BhattDist = zeros(N);
Out.KLDivergence = zeros(N);
Out.CellOverSeqRatio = zeros(N);
Out.MorisitaHorn = zeros(N);
for r = 2:N
    for c = 1:r-1
        Out.SeqOverlap(r, c) = numel(vertcat(Lookup{r, c})) / (numel(Data{r}) + numel(Data{c}) - numel(vertcat(Lookup{r, c}))) * 100; %In percentage 
        Out.SeqOverlap(c, r) = Out.SeqOverlap(r, c); 
        
        WeightSeqOverlapCt = min([vertcat(LookupWeight{r, c}), vertcat(LookupWeight{c, r})], [], 2); %Need min, as you're doing overlap of cells
        Out.CellOverlap(r, c) = sum(WeightSeqOverlapCt) / (sum(Weight{r}) + sum(Weight{c}) - sum(WeightSeqOverlapCt)) * 100; %In percentage
        Out.CellOverlap(c, r) = Out.CellOverlap(r, c);
        
        Out.BhattDist(r, c) = calcBhattStat(LookupWeight{r, c}, LookupWeight{c, r});
        Out.BhattDist(c, r) = Out.BhattDist(r, c);
        [Out.KLDivergence(r, c), Out.KLDivergence(c,r)] = calcKullbackLeibler(LookupWeight{r, c}, LookupWeight{c, r});
        
        Out.MorisitaHorn(r, c) = calcMorisitaHorn(LookupWeight{r, c}, LookupWeight{c, r});
        Out.MorisitaHorn(c, r) = Out.MorisitaHorn(r, c);
    end
end
Out.CellOverSeqRatio = Out.CellOverlap ./ Out.SeqOverlap;
