%mapDotColor_CDR3 will assign a RGB color each CDR3seq based on the
%input RefCDR3 order. Used mainly for plotTreeData.m.
%
%  [DotColor, UnqDotColor] = mapDotColor_CDR3(CDR3seq,UnqCDR3seq)
%
%  [DotColor, UnqDotColor] = mapDotColor_CDR3(TreeData,TreeHeader,UnqCDR3seq)
%
%  [DotColor, UnqDotColor] = mapDotColor_CDR3(...,'ColorMap',ColorMap)
%
%  INPUT
%    CDR3seq: Nx1 cell of CDR3 sequences to be labeled. CDR3seq that is not
%      in UnqCDR3seq will return a DotColor of [0 0 0].
%    UnqCDR3seq: Mx1 cell containing unique CDR3 sequences, ordered the
%      same way as CDR3legend output from makeTreeLegend_CDR3.m. If user
%      explicity uses UnqCDR3seq = [], then this will determine the default
%      UnqCDR3seq returned by makeTreeLegend_CDR3.m.
%    ColorMap: Qx3 matrix containing the RGB colormap that the user wants
%      to use. Default is ColorMap = jet;
%
%  OUTPUT
%    DotColor: Nx3 matrix storing RGB data for each of CDR3seq
%    UnqDotColor: Mx3 matrix strong RGB data for each of UnqCDR3seq


function [DotColor, UnqDotColor] = mapDotColor_CDR3(varargin)
%Setup the user-specified or default colormap to use
ColorMapLoc = findCell(varargin,'ColorMap');
if max(ColorMapLoc) > 0
    ColorMap = varargin{ColorMapLoc+1};
    if size(ColorMap,2) ~= 3
        error('MapDotColor_CDR3: ColorMap should be an RGB matrix');
    end
    varargin(ColorMapLoc:ColorMapLoc+1) = []; %Remove to allow for proper parsing below
else %Use the default jet colormap, mut make middle area darker (hard to see otherwise).
    ColorMap = jet;
    MidRange1 = round(size(ColorMap,1)/4);
    MidRange2 = MidRange1*3;
    ColorMap(MidRange1:MidRange2,:) = ColorMap(MidRange1:MidRange2,:)*0.6; %Making mid range color darker
end

%Determine other input types
if length(varargin) == 2
    CDR3seq = varargin{1};
    UnqCDR3seq = varargin{2};
elseif length(varargin) == 3
    TreeData = varargin{1};
    TreeHeader = varargin{2};
    TH = getTreeHeaderVar(TreeHeader);
    CDR3seq = TreeData(:,TH.CDR3seqLoc);
    UnqCDR3seq = varargin{3};
else
    error('mapDotColor_CDR3: Number of inputs is not correct');
end

%Make sure lengths are all the same
CDR3length = length(CDR3seq{1});
for j = 2:length(CDR3seq)
    if length(CDR3seq{j}) ~= CDR3length
        error('mapDotColor_CDR3: Require same-length CDR3 sequences');
    end
end

%Find the default UnqCDR3seq if user inputted an empty bracket UnqCDR3seq
if isempty(UnqCDR3seq)
    [~,UnqCDR3seq,~] = makeTreeLegend_CDR3(CDR3seq);
end
    
%Get the unique color for reference sequence, if not provided    
ColorIdx = round(linspace(1,size(ColorMap,1),length(UnqCDR3seq)));
if length(UnqCDR3seq) == 1
    ColorIdx = 1; %Always make it first color for singletons
end
UnqDotColor = ColorMap(ColorIdx,:);

%Identify which CDR3seq belongs to what RefCDR3.
CDR3map = zeros(size(CDR3seq,1),1);
for j = 1:size(CDR3seq,1)
    for k = 1:size(UnqCDR3seq,1)
        if strcmpi(CDR3seq{j},UnqCDR3seq{k})
            CDR3map(j) = k;
            break
        end
    end
end

%Assign a Dot color for every CDR3seq
DotColor = zeros(size(CDR3seq,1),3);
DotColor(CDR3map>0,:) = UnqDotColor(CDR3map(CDR3map>0),:);


