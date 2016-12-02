%combineImages will take a list of image files, and then concantenate them.
%Useful for making collages. Will align by top.

function ImageC = combineImages(varargin)
%See if "SaveName" is provided
SaveName = '';
length(varargin)
for j = 1:length(varargin)-1
    if strcmpi(varargin{j},'savename')
        SaveName = varargin{j+1};
        break;
    end
end
varargin(j:j+1) = [];

if isempty(varargin)
    [FileNames, FilePath] = uigetfile('*.png','Select image files','multiselect','on');
    if ischar(FileNames)
        error('Only works for multiple image file')
    end
else
    FileNames = varargin{1};
    if nargin >= 2
        FilePath = varargin{2};
    else
        if ispc
            FilePath = [cd '\'];
        else
            FilePath = [cd '/'];
        end
    end
end

AllIm = cell(1,length(FileNames));
MaxH = 0;
MaxW = 0;
MaxZ = 0;
for f = 1:length(FileNames)
    CurIm = imread([FilePath FileNames{f}]);
    CurH = size(CurIm,1);
    CurW = size(CurIm,2);
    CurZ = size(CurIm,3);
    
    %Determine max H
    if CurH > MaxH
        MaxH = CurH;
    end
    
    %Determine max H
    if CurZ > MaxZ
        MaxZ = CurZ;
    end
    
    %Continuously add W;
    MaxW = MaxW + CurW;
    
    AllIm{f} = CurIm;
end

%Create the concantenated image files
ImClass = class(AllIm{1});
ImageC = zeros(MaxH,MaxW,MaxZ,ImClass);
W1 = 1; %Keep track of W start position.
for f = 1:length(FileNames)
    CurIm = AllIm{f};
    CurH = size(CurIm,1);
    CurW = size(CurIm,2);
    CurZ = size(CurIm,3);
    
    W2 = W1 + CurW - 1;
    ImageC(1:CurH,W1:W2,1:CurZ) = CurIm;
    W1 = W2 + 1;
end

if isempty(SaveName)
    [SaveName, SavePath] = uiputfile('*.png','Save combined image as');
    imwrite(ImageC,[SavePath SaveName]);
else
    imwrite(ImageC,[SaveName]);
end

