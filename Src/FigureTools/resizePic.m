%resizePic will take an image file and resize it to a different size and
%format.
%
%  resizePic
%
function resizePic
[FileName, FilePath] = uigetfile('*.pdf;*.png;*.tif;*.jpg', 'Open image file');
DotIdx = find(FileName == '.', 1, 'last');
FilePre = FileName(1:DotIdx-1);
FileExt = FileName(DotIdx:end);

Gx = figure;
Ax = axes;
Im = imread(fullfile(FilePath, FileName));
imshow(Im, 'parent', Ax);

OKAY = 0;
while ~OKAY
    DPI = round(input('What is the desired DPI? '));
    W = input('What is the desired width in inches? ');
    H = input('What is the desired height in inches? ');

    resizeFigure(Gx, W, H)
    Ax.Units = 'normalized';
    Ax.Position = [0 0 1 1];
    
    OKAY = isempty(input('Enter to continue, or ''n'' to redo: ', 's'));
end

SaveName = fullfile(FilePath, sprintf('%s_H%0.2f_W%0.2f_DPI%d%s', FilePre, H, W, DPI, FileExt));
[SaveName, SavePath] = uiputfile('*.*', 'Save image as', fullfile(FilePath, SaveName));
savePlot(Gx', 'saveas', fullfile(SavePath, SaveName), 'DPI', DPI)
close(Gx);

