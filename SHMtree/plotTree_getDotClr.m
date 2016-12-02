%plotTree_getDotClr will assign a RGB color each CDR3seq, based on the
%input RefCDR3. DotClr is the RGB for each CDR3seq, and UnqClr is the color
%fo each unique CDR3seq. We need UnqClr for making colored text legend.
%
%  [DotClr,UnqClr] = plotTree_getDotClr(CDR3seq,RefCDR3)   will assign a jet color
%  to each CDR3seq, based on the order of RefCDR3. Use
%  plotTree_getUnqCDR3.m before to get RefCDR3 (or UnqCDR3 output), sorted
%  by ham distance to root seq.
%
%  DotClr = plotTree_getDotClr(CDR3seq,RefCDR3,RefClr)   will use int
%  RefClr fo each RefCDR3, and assign this to the CDR3seq instead.
%
%  

function [DotClr, RefClr] = plotTree_getDotClr(CDR3seq,varargin)
P = inputParser;
addRequired(P,'CDR3seq',@iscell);
addOptional(P,'RefCDR3',[],@iscell);
addOptional(P,'RefClr',[],@isnumeric);
parse(P,CDR3seq,varargin{:});

CDR3seq = P.Results.CDR3seq;
RefCDR3 = P.Results.RefCDR3;
RefClr = P.Results.RefClr;

%==========================================================================
%Determine the color scheme, but use custom color if it exists

if isempty(RefCDR3)
    [RefCDR3,~,~] = plotTree_getUnqCDR3(CDR3seq);
end
    
%Get the unique color for reference sequence, if not provided
if isempty(RefClr)
    ClrMap = jet;
    MidRange1 = round(size(ClrMap,1)/4);
    MidRange2 = MidRange1*3; 
    ClrMap(MidRange1:MidRange2,:) = ClrMap(MidRange1:MidRange2,:)*0.6; %Making mid range color darker
    
    ClrIdx = round(linspace(1,size(ClrMap,1),length(RefCDR3)));
    if length(RefCDR3) == 1
        ClrIdx = 1; %Always make it first color for singletons
    end
    RefClr = ClrMap(ClrIdx,:);
end

%Identify which CDR3seq belongs to what RefCDR3.
CDR3map = zeros(size(CDR3seq,1),1);
for j = 1:size(CDR3seq,1)
    for k = 1:size(RefCDR3,1)
        if strcmpi(CDR3seq{j},RefCDR3{k})
            CDR3map(j) = k;
            break
        end
    end
end

%Assign a Dot color for every CDR3seq
DotClr = zeros(size(CDR3seq,1),3);
DotClr(CDR3map>0,:) = RefClr(CDR3map(CDR3map>0),:);


