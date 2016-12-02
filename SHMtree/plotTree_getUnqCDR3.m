%plotTree_getUnqCDR3 will sort the CDR3seq cell string, return the unique
%CDR3 based on  hamming distance from the root CDR3. If AncMap (ancestral
%map matrix) is provided, it will determine the root using AncMap. The
%output is a sorted cell string for unique CDR3s.
%
%  [UnqCDR3, UnqIdx, UnqCDR3align] = plotTree_getUnqCDR3(CDR3seq)   will
%  assume first CDR3seq is root.
%
%  [UnqCDR3, UnqIdx, UnqCDR3align] = plotTree_getUnqCDR3(CDR3seq,AncMap)
%  will use AncMap to find root CDR3.
%
%  UnqCDR3 = unique CDR3 of CDR3seq, ordered from root to furtherest from
%  root.
%
%  UnqIdx = matrix of number that corresponds each CDR3seq to the UnqSeq
%  index. This is the same size as CDR3seq.
%
%  UnqCDR3align is a character array (not a cell) with all CDR3 alignment,
%  replacing matched letters with '-'. This is used later to build plotTree
%  legends.

function [UnqCDR3, UnqIdx, UnqCDR3align] = plotTree_getUnqCDR3(CDR3seq,varargin)
P = inputParser;
addRequired(P,'CDR3seq',@iscell);
addOptional(P,'AncMap',[],@isnumeric);
parse(P,CDR3seq,varargin{:});

CDR3seq = P.Results.CDR3seq;
AncMap = P.Results.AncMap;

%Determine unique CDR3s, sorted by distance to root CDR3
if ~isempty(AncMap)   
    RootLoc = find(AncMap(:,2) == 0);
    RootLoc = RootLoc(1);
else
    RootLoc = 1;
end
RootCDR3 = CDR3seq{RootLoc};

%Finding unique ones and then determining hamming distance to RootCDR3
[UnqCDR3,~,UnqIdx] = unique(CDR3seq);
SortMat = zeros(size(UnqCDR3,1),1);
for k = 1:length(UnqCDR3)
    SortMat(k,1) = sum(UnqCDR3{k} == RootCDR3);   
end
[~,SortIDX] = sort(SortMat,'descend');
UnqCDR3 = UnqCDR3(SortIDX); %Output1

%Determine the new numbering for UnqIdx, now that you resorted UnqCDR3, and
%append this to AncMap
UnqIdxT = UnqIdx;
for k = 1:length(UnqCDR3)
    UnqIdx(UnqIdxT == SortIDX(k)) = k;
end

%Determine the legend text here, simplifying matches by '_'
UnqCDR3align = char(UnqCDR3);
MatchLoc = zeros(size(UnqCDR3align))>1;
for j = 2:size(UnqCDR3align,1)
    MatchLoc(j,:) = RootCDR3 == UnqCDR3align(j,:);
end
UnqCDR3align(MatchLoc) = '-';