% Out
%   .KLDivergence %Kullback_Leibler
%   .BhattDist    %Bhattacharya Distance matrix
%   .Lookup  % creates a NxN cell of Row vs Col repertoire Idx of overlap. This can further be grouped according to filters.
%     .Clone %Idx of clones with same CDR3 / total unique clones (no template data used)
%     .AllClonotype %Idx of clones with same CDR3 AND single-member clonotypes / total clonotypes
%     .SingleClonotype %Idx of clones with same CDR3 AND single-member clonotypes / total clonotypes
%     .BranchedClonotype %Idx of clones with same CDR3 AND NOT single-member clonotypes / total clonotypes
%     .Cell %Idx of cells with same CDR3 / total number of cells (template data used)
% 
% 
% 
% 
% 
% 
% 
% 

function Heatmap = getConvData(Data)
[UnqSeq, UnqIdx] = getUnqSeq(Data);
Lookup = getLookup(Data, UnqSeq, UnqIdx);
Template = getTemplate(Data, Lookup);
GrpNum = getGrpNum(Data, Lookup);
Heatmap = getHeatmap(Data, Lookup, UnqSeq, UnqIdx, Template, GrpNum);

function [UnqSeq, UnqIdx] = getUnqSeq(Data)
N = numel(Data);
UnqSeq = cell(N, 1);
UnqIdx = cell(N, 1);
for f = 1:N
    [UnqSeq{f}, ~, ~, UnqIdx{f}] = unique2(Data(f).hCDR3);
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

function Template = getTemplate(Data, Lookup)
N = numel(Data);
Template = cell(N);
for r = 2:N
    TemplateR = Data(r).Template;
    for c = 1:r-1
        TemplateC = Data(c).Template;
        Template{r, c} = cellfun(@(x) TemplateR(x), Lookup{r, c}, 'un', 0);
        Template{c, r} = cellfun(@(x) TemplateC(x), Lookup{c, r}, 'un', 0);
    end
end

function GrpNum = getGrpNum(Data, Lookup)
N = numel(Data);
GrpNum = cell(N);
for r = 2:N
    for c = 1:r-1
        GrpNum{r, c} = cellfun(@(x) Data(r).GrpNum(x), Lookup{r, c}, 'un', 0);
        GrpNum{c, r} = cellfun(@(x) Data(c).GrpNum(x), Lookup{c, r}, 'un', 0);
    end
end

function Heatmap = getHeatmap(Data, Lookup, UnqSeq, UnqIdx, Template, GrpNum)
N = numel(Data);
CloneOvlp = zeros(N);
for r = 2:N
    for c = 1:r-1
        %Clone Level = # shared CDR3 for all clone / total # unique clone (no template)
        CloneOvlp(r, c) = numel(vertcat(Lookup{r, c})) / numel(Data(r).GrpNum); 
        CloneOvlp(c, r) = numel(vertcat(Lookup{c, r})) / numel(Data(c).GrpNum); 
    end
end
        
CDR3Ovlp = zeros(N);
for r = 2:N
    for c = 1:r-1
        %CDR3 Level = # shared unique CDR3 / total # unique CDR3
        CDR3Ovlp(r, c) = numel(Lookup{r, c}) / numel(UnqSeq{r});
        CDR3Ovlp(c, r) = numel(Lookup{c, r}) / numel(UnqSeq{c});
    end
end

BhattDist = zeros(N);
KLDivergence = zeros(N);
for r = 2:N
    for c = 1:r-1
        DistR = cellfun(@(x) sum(Data(r).Template(x)), Lookup{r, c});
        DistC = cellfun(@(x) sum(Data(c).Template(x)), Lookup{c, r});

        BhattDist(r, c) = calcBhattStat(DistR, DistC);
        BhattDist(c, r) = BhattDist(r, c);

        [KLDivergence(r, c), KLDivergence(c,r)] = calcKullbackLeibler(DistR, DistC);
    end
end

Heatmap.CloneOvlp = CloneOvlp;
Heatmap.CDR3Ovlp  = CDR3Ovlp; 
Heatmap.BhattDist = BhattDist; 
Heatmap.KLDivergence = KLDivergence;


