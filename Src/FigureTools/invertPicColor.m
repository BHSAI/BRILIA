%invertPicColor will invert image colors of selected files in the computer.
%
%  invertPicColor
%
%  invertPicColor(FileNames)
%
%  INPUT
%    FileNames: full file name of an image file to invert the color of. 
%      If empty, will ask use to select files. 
%      If char, will look for file(s) such as C:\Images\*.png.
%      If cell, will look for file(s) in cell array as is. 
%
%  OUTPUT
%    Saves the inverted color pictures in the same directory as the
%    original file, but with a ".inv" added before the file extension. 
%
%  EXAMPLE
%    invertPicColor('C:\Test\*.png'); %Will invert al .png files in C:\Test
%                                     %folder as *.inv.png files.
%
function invertPicColor(FileNames)
if nargin == 0 || isempty(FileNames)
    [FileName, FilePath] = uigetfile('*.png;*.tif;*.jpg', 'Open picture files', 'MultiSelect', 'on');
    if isnumeric(FileName); return; end
    FileNames = fullfile(FilePath, FileName);
end
if ischar(FileNames)
    FileNames = dir(FileNames);
    FileNames = fullfile({FileNames.folder}, {FileNames.name});
end
if isempty(FileNames); return; end

for f = 1:length(FileNames)
    if ~exist(FileNames{f}, 'file')
        warning('%s: Could not find file "%s". Skipping.', mfilename, FileNames{f});
    end
    
    Im = imread(FileNames{f});
    if ismember(class(Im), {'double', 'logical'})
        MaxVal = 1;
    else
        MaxVal = intmax(class(Im));
    end
    Im = MaxVal - Im;
    
    [Path, Name, Ext] = fileparts(FileNames{f});
    SaveName = fullfile(Path, [Name '.inv' Ext]);
    imwrite(Im, SaveName)
end

fprintf('%s: Inverted %d images.\n', mfilename, length(FileNames));