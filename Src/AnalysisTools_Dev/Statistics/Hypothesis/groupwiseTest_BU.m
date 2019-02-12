%groupwiseTest is a procedural test looking for difference between
%distributions of categorical variables via bootstrapping difference in
%observed values, pairwise comparison of each bar frequencies, groupwise
%comparison of observed differences, multiple hypothesis correction, and
%p-value calculation of probability that an observed difference in groups
%is caused by random error alone. 
%
%  [Pmat, Gmat, Mmat] = groupwiseTest(Data, Group, varargin)
%
%  INPUT
%    Data: MxN vector of frequencies (>= 0 value) for M observations and N
%      distributions.
%    Group: 1xN vector of integers assigning each data column to a group
%    Param      Value      Description
%    ---------- ---------- ------------------------------------------------
%    Size       N > 0      The number of samples to draw with replacement 
%                          from each nth column distribution. Default is
%                          min(sum(Data, 1)), or the smallest count for all
%                          distributions.
%    Iter       N > 0      The number of resampling iterations to perform.
%                          Recommend > 100 to get 95% confidence intervals.
%                          Default: 100
%    Alpha      0<Alpha<1  The significance level for the conf. interval.
%
%  OUTPUT
%    Sall: M(Nx(N+2) vector of P values from pairwiseTest for M variables,
%      where the first 2 col are Sub1# and Sub2# for pairwise comparison.


function [Pall, Gall, MED, CI, ScaleMat] = groupwiseTest_BU(Data, Group, varargin)
P = inputParser;
P.KeepUnmatched = 1;
addParameter(P, 'Alpha', 0.05, @(x) x >= 0 && x < 1);
parse(P, varargin{:});
Alpha = P.Results.Alpha;

if iscell(Data)
    Data = cell2mat(Data);
end
[M, N] = size(Data);

%Bootstrap to get the CI attributed by sampling error
[MED, CI, ~, ~, BootData] = bootstrapDist(Data, varargin{:});

%Rescale the frequencies because the bootstrap would have affected the
%total population of a frequency and thus needs to be readjusted
ScaleMat = zeros(1, N);
for j = 1:N
    [MED(:, j), ScaleMat(j)] = rescaleDist(Data(:, j), MED(:, j));
    BootData{j} = ScaleMat(j) * BootData{j};
    CI(:, j) = cellfun(@(x) ScaleMat(j) * x, CI(:, j), 'un', 0);
end

%Compute groupwise differences
[UnqGroup, ~, UnqIdx] = unique(Group);
G = length(UnqGroup);
UnqCellIdx = cell(1, length(UnqGroup));
for y = 1:length(UnqGroup)
    UnqCellIdx{y} = find(UnqIdx == y);
end

Pvec = cell(1, M);
Gvec = cell(1, M);
CombVec = combvec2(1:size(BootData{1}, 2), 1:size(BootData{1}, 2));
for m = 1:M
    %Compute pairwise differences
    Pmat = zeros(N);
    for s1 = 1:N
        for s2 = 1:s1-1
            [~, Pmat(s1, s2)] = pairwiseTest(BootData{s1}(m, :), BootData{s2}(m, :), 'Alpha', Alpha, 'CombVec', CombVec);
            Pmat(s2, s1) = Pmat(s1, s2);
        end
    end
    
    %Inter group level comparison of significance, using BH correction.
    Gmat = zeros(G);
    for g1 = 1:G
        for g2 = 1:g1-1
            Gall = Pmat(UnqCellIdx{g1}, UnqCellIdx{g2});
            SigIdx = correctP(Gall(:), Alpha, 'bh');
            H = zeros(size(Gall));
            H(SigIdx) = 1;
            Gmat(g1, g2) = binormCDF(sum(H(:))-1, numel(H), 0.5, 'upper');
            Gmat(g2, g1) = Gmat(g1, g2);
        end
    end
    
    %Format results in printable format
    Pvec{m} = squareform(Pmat, 'tovector')';
    Gvec{m} = squareform(Gmat, 'tovector')';
end

Pall = [nchoosek(1:N, 2) Pvec{:}]; %Sub1# Sub2# Var1 Var2 ... VarN | 1 (sig diff) or 0 (no diff)
Gall = [nchoosek(1:G, 2) Gvec{:}]; %Sub1# Sub2# Var1 Var2 ... VarN | binocdf(#Sig, #Possible, 0.5, 'upper')
