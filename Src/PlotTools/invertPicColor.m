%invertPicColor will invert image colors but will. This is used mainly if
%you want to show white-background images on a black-background powerpoint,
%or vice versa.
%
%  invertPicColor
%
%  invertPicColor(FullFileName)
%
%  INPUT
%    FullFileName: file path and name of image file to invert the color. If
%      left empty, will ask use to seleect files. FullFileName can be
%      either a char or cell of char file names.
%
%  NOTE
%    Automatically saves the inverted pictures in the same directory as
%    the picture files, with a ".inv" added before the file extension.

function invertPicColor(varargin)
if isempty(varargin)
    [FileName, FilePath] = uigetfile('*.png;*.tif;*.jpg','Open picture files','MultiSelect','on');
else
    FileName = varargin{1};
end
if ischar(FileName)
    FileName = {FileName};
end

if ~isempty(regexp(FileName{1}, filesep, 'once'))
   FullFileName = FileName;
else
    if ~exist('FilePath','var')
        FilePath = [cd filesep];
    end
    FullFileName = cell(length(FileName),1);
    for f = 1:length(FileName)
        FullFileName{f} = [FilePath FileName{f}];
    end
end

for f = 1:length(FullFileName)
    [FilePath, FileName, FileExt] = parseFileNameT(FullFileName{f});
    
    Im = imread(FullFileName{f});
    switch class(Im)
        case 'uint8'
            MaxVal = 2^8-1;
        case 'uint16'
            MaxVal = 2^16-1;
        case 'uint32'
            MaxVal = 2^32-1;
        case 'double'
            MaxVal = 1;
        case 'logical'
            MaxVal = 1;
    end
    
    Im = MaxVal - Im;
    DotLoc = find(FileName == '.');
    SaveName = [FileName(1:DotLoc(end)) 'inv' FileExt];
    imwrite(Im, [FilePath SaveName])
end




