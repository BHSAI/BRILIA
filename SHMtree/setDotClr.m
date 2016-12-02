%setDotClr will assign a dot color associated each CDR3 seq.  Otherwise, It will use the given dot color, and use Jet color
%schemes to find equally-spaced colors for each unique CDR3. If no
%UnqCDR3seq is provided, then this assumes CDR3seq(1) is the "root", and
%will according define colors such that root is blue, and most distance
%away from root is red.
%  setDotClr(CDR3seq,UnqCDR3seq,UnqClr)   If reference UnqCDR3 and
%  UnqCDR3Color are provided, it will match the dot color accordingly.

function [DotClr, vararout] = setDotClr(CDR3seq,varargin)
P = inputParser;
addRequired(P,'CDR3seq',@iscell);
addOptional(P,'UnqCDR3seq','',@ischar);
addOptional(P,'UnqClr',[],@isnumeric);
parseInput(P,CDR3seq,varargin{:});

CDR3seq = P.results.CDR3seq;
UnqCDR3seq = P.results.UnqCDR3seq;
UnqClr = P.results.UnqClr;

%==========================================================================
%Determine unique CDR3s, sorted by distance to root CDR3
RootLoc = find(AncMap(:,2) == 0);
RootLoc = RootLoc(1);
RootCDR3 = CDR3seq{RootLoc};
[UnqCDR3,~,UnqIdx] = unique(CDR3seq);
SortMat = zeros(size(UnqCDR3,1),1);
for k = 1:length(UnqCDR3)
    SortMat(k,1) = sum(UnqCDR3{k} == RootCDR3);   
end
[~,SortIDX] = sort(SortMat,'descend');
UnqCDR3 = UnqCDR3(SortIDX(:,1));

%Determine the new numbering for UnqIdx, now that you resorted UnqCDR3, and
%append this to AncMap
UnqIdxT = UnqIdx;
for k = 1:length(UnqCDR3)
    UnqIdx(UnqIdxT == SortIDX(k)) = k;
end
AncMap = [AncMap UnqIdx];

    %Determine the legend text here, simplifying matches by '_'
    AllCDR3 = char(UnqCDR3);
    MatchLoc = zeros(size(AllCDR3))>1;
    for j = 2:size(AllCDR3,1)
        MatchLoc(j,:) = RootCDR3 == AllCDR3(j,:);
    end
    AllCDR3(MatchLoc) = '-';
else
    UnqCDR3 = cell(size(AncMap,1),1);
end









%Determine the color scheme, but use custom color if it exists
if isempty(varargin)
    %Set the unique CDR color map, but darken some that are too bright
    ClrMap = colormap('jet');
    MidRange1 = round(size(ClrMap,1)/4);
    MidRange2 = MidRange1*3; 
    ClrMap(MidRange1:MidRange2,:) = ClrMap(MidRange1:MidRange2,:)*0.6;

    %Give CDR3s unique color codes
    if length(UnqCDR3) > 1
        ClrIdx = round(linspace(1,size(ClrMap,1),length(UnqCDR3)));
    else
        ClrIdx = 1; %Always make it first color for singletons
    end
    UnqClr = ClrMap(ClrIdx,:);

    %Assign a DotColor for every entity in AncMap
    DotClr = zeros(size(AncMap,1),3);
    for j = 1:size(AncMap,1)
        DotClr(j,:) = UnqClr(AncMap(j,end),:);
    end
end


