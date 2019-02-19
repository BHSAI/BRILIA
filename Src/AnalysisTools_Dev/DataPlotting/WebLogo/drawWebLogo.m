%drawWebLogo takes an aligned amino acid sequence and returns a web logo
%picture.
%   
%  IMlogo = drawWebLogo(Seq)
%
%  IMlogo = drawWebLogo(Seq, Range)
%
%  [IMlogo, IMmap] = drawWebLogo(Seq,...) 
%
%  INPUT
%    Seq: set of an amino acid sequences to draw a web logo for
%    Range: range of sequences to draw the web logos for 
%
%  OUTPUT
%    IMlogo: web logo image
%    IMmap: matrix of size IM indicating the Nth position of the sequence.
%      This is used mainly to map the sequence position to the position the
%      user clicks on the web logo image IM.
%
%  Example:
%        Seq = {'LSGGQRQRVAIARALAL'; 
%               'LSGGEKQRVAIARALMN'; 
%               'LSGGQIQRVLLARALAA';
%               'LSGGERRRLEIACVLAL'; 
%               'FSGGEKKKNELWQMLAL'; 
%               'LSGGERRRLEIACVLAL'};
%        IMlogo = drawWebLogo(Seq,1:3);
%        imshow(IMlogo)

function varargout = drawWebLogo(Seq, varargin)
if iscell(Seq)
    Seq = cell2mat(Seq);
end
if numel(varargin) == 1
    Range = varargin{1};
else
    Range = 1:size(Seq, 2);
end

FontPath = fullfile(fileparts(mfilename('fullfile')), 'Font');

%Generate weblogo
LogoH = 600; %Pixel height of logo
FontW = 300; %Pixel width of each letter

IMlogo = [];
for kk = 1:length(Range)

    PosNum = Range(kk);
    AaFreq = cell2mat(struct2cell(aacount(Seq(:, PosNum))));
    KeepLoc = (AaFreq > 0);
    AaKeep = int2aa(find(KeepLoc == 1));
    AaFreqNorm = AaFreq(KeepLoc)/sum(AaFreq);
    
    %Sort, largest item top
    AaAll = sortrows([AaFreqNorm (1:length(AaKeep))']);
    AaFreqNorm = AaAll(:,1);
    AaKeep = AaKeep(AaAll(:,2));

    IMcol = [];
    FontHeightTot = 0;
    for qq = 1:length(AaKeep)
        FontFile = fullfile(FontPath, [AaKeep(qq) '.png']);
        IM1 = imread(FontFile);
        FontHeight = round(LogoH*AaFreqNorm(qq));
        if qq == length(AaKeep)
            FontHeight = LogoH - FontHeightTot;
        end
        if FontHeight < 1
            continue
        end
        FontHeightTot = FontHeightTot+FontHeight;
        IM1 = shrinkImage(im2double(IM1), FontHeight, FontW);
        if isempty(IMcol)
            IMcol = IM1;
        else
            IMcol = cat(1,IMcol,IM1);
        end
    end
    if isempty(IMlogo)
        IMlogo = IMcol;
    else
        IMlogo = cat(2,IMlogo,IMcol);
    end
end

%Generate clickermap
PosLen = length(Range);
Width = size(IMlogo,2)/PosLen;
IMmap = zeros(size(IMlogo));
for j = 1:PosLen
    IMmap(:,1+(j-1)*Width:j*Width) = j;
end

varargout{1} = IMlogo;
varargout{2} = IMmap;