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
%    Size: the number of samples to draw with replacement from each nth
%      column distribution. Default: min(sum(Data, 1))
%    Iter: the number of iterations to perform. Default: 100
%    Alpha: the significance level for the confidence interval
%
%  OUTPUT
%    Sall: M(Nx(N+2) vector of P values from pairwiseTest for M variables,
%      where the first 2 col are Sub1# and Sub2# for pairwise comparison.

function [Pall, Gall, MED, CI, ScaleMat] = groupwiseTest(Data, Group, varargin)
P = inputParser;
P.KeepUnmatched = 1;
addParameter(P, 'Alpha', 0.05, @(x) x >= 0 && x < 1);
addParameter(P, 'Prefix', '' , @ischar);
parse(P, varargin{:});
Alpha = P.Results.Alpha;
Prefix = P.Results.Prefix;

%Temp code for showing how this works
PrintThis = 0;
PrintIdx = findCell(varargin, 'PrintThis');
if PrintIdx > 1
    PrintThis = 1;
    VarName = varargin{PrintIdx+1};
    varargin(PrintIdx:PrintIdx+1) = [];
else
    VarName = [];
end
%End of Temp code

if iscell(Data)
    Data = cell2mat(Data);
end
[M, N] = size(Data);

%Perform boot strap to get the 95% CI attributed by sampling error
[MED, CI, ~, ~, BootData] = bootstrapDist(Data, varargin{:});

%Determine the scaling matrix for sampling population control
ScaleMat = zeros(1, N);
for j = 1:N
    [MED(:, j), ScaleMat(j)] = rescaleDist(Data(:, j), MED(:, j));
    for k = 1:M
        CI{k, j} = ScaleMat(j) * CI{k, j};
    end
    BootData{j} = ScaleMat(j) * BootData{j};
end

%For each cateogry, compute pairwise differences
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
    %Subject level pairwise comparison
    Pmat = zeros(N);
    for s1 = 1:N
        for s2 = 1:s1-1
            [~, Pmat(s1, s2)] = pairwiseTest(BootData{s1}(m, :), BootData{s2}(m, :), 'CombVec', CombVec, 'Alpha', Alpha);
            if Pmat(s1, s2) > 1; pause; end
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
%     
%     if ~isempty(VarName) && contains(VarName{m}, {'IGHV8-8*01', 'IGHV1-64*01'}) && PrintThis
%         warning('%s: Don''t forget about removing this msg!', mfilename);
%         
%         PrintCI = vertcat(CI{m, :})';
%         
%         plotErrorBox(MED(m, :), 'CI', CI(m, :));
%         set(gca, 'XTick', 1:5, 'XTickLabel', VarName(m));
%         resizeFigure(gcf, 3, 6);
%         drawnow
%         resizeSubplots;
%         
%         SaveNamePre = strrep([Prefix '_' VarName{m}], '*', '-');
%         savePlot(gcf, 'SaveAs', ['Example_' SaveNamePre '.png']);
%         
%         writeDlmFile(num2cell([1:size(Pmat, 1); Group; PrintCI(1, :); MED(m, :); PrintCI(2, :)]), [SaveNamePre '_MED.csv'], ',');
%         writeDlmFile(num2cell([1:size(Pmat, 1); Group; Pmat]), [SaveNamePre '_Pmat.csv'], ',');
%         writeDlmFile(num2cell([UnqGroup; Gmat]), [SaveNamePre '_Gmat.csv'], ',');
%     end
    
    %Format results in printable format
    Pvec{m} = squareform(Pmat, 'tovector')';
    Gvec{m} = squareform(Gmat, 'tovector')';
end

Pall = [nchoosek(1:N, 2) Pvec{:}]; %Sub1# Sub2# Var1 Var2 ... VarN | 1 (sig diff) or 0 (no diff)
Gall = [nchoosek(1:G, 2) Gvec{:}]; %Sub1# Sub2# Var1 Var2 ... VarN | binocdf(#Sig, #Possible, 0.5, 'upper')
