%invertPicColor will invert image colors of selected files in the computer.
%
%  invertPicColor
%
%  invertPicColor(FileNames)
%
%  INPUT
%    FileNames: full file name of an image file to invert the color of. If
%      left empty, will ask use to select files. FileNames can be
%      either a char (ex: 'image1.jpg' or a cell array of file names 
%      (ex: {'image1.jpg' 'image2.png'}).
%
%  OUTPUT
%    Automatically saves the inverted pictures in the same directory as
%      the picture files, with a ".inv" added before the file extension.

function invertPicColor(FileNames)
%Image file(s) selection
if nargin == 0 || isempty(FileNames)
    [FileName, FilePath] = uigetfile('*.png;*.tif;*.jpg', 'Open picture files', 'MultiSelect', 'on');
    if isnumeric(FileName) 
        return 
    elseif ischar(FileName)
        FileName = {FileName};
    end
    for f = 1:length(FileName)
        FileNames{f} = fullfile(FilePath, FileName{f});
    end
elseif ischar(FileNames)
    FileNames = {FileNames};
end

%Image inversion
for f = 1:length(FileNames)
    if ~exist(FileNames{f}, 'file')
        error('%s: Could not find file "%s".', mfilename, FileNames{f});
    end
    
    Im = imread(FileNames{f});
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
    
    [Path, Name, Ext] = fileparts(FileNames{f});
    SaveName = fullfile(Path, [Name '.inv' Ext]);
    imwrite(Im, SaveName)
end

fprintf('%s: Inverted %d images.\n', mfilename, length(FileNames));