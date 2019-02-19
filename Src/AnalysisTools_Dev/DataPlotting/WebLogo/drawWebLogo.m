%drawWebLogo takes an aligned amino acid sequence and returns a web logo
%picture.
%   
%  IMlogo = drawWebLogo(Seq)
%
%  IMlogo = drawWebLogo(Seq,PosNums)
%
%  [IMlogo, IMmap] = drawWebLogo(Seq,...) 
%
%  INPUT
%    Seq: set of an amino acid sequences to draw a web logo for
%    PosNums: range of sequences to draw the web logos for 
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

function varargout = drawWebLogo(Seq,varargin)
if iscell(Seq)
    Seq = cell2mat(Seq);
end
if length(varargin) == 1
    PosNums = varargin{1};
else
    PosNums = 1:size(Seq,2);
end

FontPath = fullfile(fileparts(mfilename('fullfile')), 'Font');

%Generate weblogo
H = 200; %pixels
Wf = 100; %Width per font.

IMlogo = [];
for kk = 1:length(PosNums)

    PosNum = PosNums(kk);
    AAFreq = cell2mat(struct2cell(aacount(Seq(:,PosNum))));
    KeepStuff = (AAFreq > 0);
    AAkeep = int2aa(find(KeepStuff == 1));
    AAfreq = AAFreq(KeepStuff)/sum(AAFreq);
    
    %Sort, largest item top
    AAall = sortrows([AAfreq (1:length(AAkeep))']);
    AAfreq = AAall(:,1);
    AAkeep = AAkeep(AAall(:,2));

    IMcol = [];
    FontHeightTot = 0;
    for qq = 1:length(AAkeep)
        FontFile = fullfile(FontPath, [AAkeep(qq) '.png']);
        IM1 = imread(FontFile);
        FontHeight = round(H*AAfreq(qq));
        if qq == length(AAkeep)
            FontHeight = H - FontHeightTot;
        end
        if FontHeight < 1
            continue
        end
        FontHeightTot = FontHeightTot+FontHeight;
        IM1 = shrinkImage(im2double(IM1), FontHeight, Wf);
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
PosLen = length(PosNums);
Width = size(IMlogo,2)/PosLen;
IMmap = zeros(size(IMlogo));
for j = 1:PosLen
    IMmap(:,1+(j-1)*Width:j*Width) = j;
end

if nargout >= 1
    varargout{1} = IMlogo;
    if nargout == 2
        varargout{2} = IMmap;
    end
end