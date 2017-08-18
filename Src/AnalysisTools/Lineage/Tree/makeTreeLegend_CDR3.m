%makeTreeLegend_CDR3 will take CDR3 sequences in a Mx1 cell, and return a
%cell matrix of unique CDR3s that have been compared against the first
%CDR3. This is used mainly to provide the legend text for trees.
%
%  [CDR3legend, UnqCDR3seq, Idx] = makeTreeLegend_CDR3(CDR3seq)
%
%  INPUT
%    CDR3seq: Mx1 cell of CDR3 amino acid sequences, where the first one is
%      considered the root. 
%
%  OUTPUT
%    CDR3legend: a Mx1 cell containing formatted CDR3 sequences, where
%      matched letters to the first CDR3seq are replaced by '-'.
%    UnqCDR3seq: a Mx1 cell containing the unique CDR3 sequences, ordered
%      the same way as CDR3legend.
%    Idx: a Nx1 matrix same length as CDR3seq that maps the Mth position of
%      UnqCDR3seq to each entry of CDR3seq. This index is used to help
%      assign a unique color to each unique CDR3 sequence.

function [CDR3legend, UnqCDR3seq, Idx] = makeTreeLegend_CDR3(varargin)
CDR3legend = [];
UnqCDR3seq = [];
Idx = [];
if isempty(varargin)
    return;
end
CDR3seq = varargin{1};

%Make sure lengths are all the same
CDR3length = length(CDR3seq{1});
for j = 2:length(CDR3seq)
    if length(CDR3seq{j}) ~= CDR3length
        error('makeTreeLegend_CDR3: Require same-length CDR3 sequences');
    end
end

%Finding unique ones and then determining hamming distance to RootCDR3
RootCDR3 = CDR3seq{1};
[UnqCDR3seq, ~, Idx] = unique(CDR3seq);
HamDistMat = zeros(size(UnqCDR3seq, 1), 1);
for k = 1:length(UnqCDR3seq)
    HamDistMat(k, 1) = sum(UnqCDR3seq{k} == RootCDR3);   
end
[~, SortIDX] = sort(HamDistMat, 'descend');
UnqCDR3seq = UnqCDR3seq(SortIDX);

%Determine the new numbering for Idx due to resorting UnqCDR3seq
IdxT = Idx;
for k = 1:length(UnqCDR3seq)
    Idx(IdxT == SortIDX(k)) = k;
end

%Determine the legend text here, replacing matched letters with '-'
CDR3legend = cell(size(UnqCDR3seq));
CDR3legend{1} = RootCDR3;
for j = 2:length(CDR3legend)
    CurCDR3seq = UnqCDR3seq{j};
    DivideLoc = CurCDR3seq == ':';
    MatchLoc = RootCDR3 == CurCDR3seq;
    CurCDR3seq(MatchLoc) = '-';
    CurCDR3seq(DivideLoc) = ':';
    CDR3legend{j} = CurCDR3seq;
end
