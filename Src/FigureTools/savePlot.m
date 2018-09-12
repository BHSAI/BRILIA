%savePlot will save the plot specified by the gcf, save name, and save
%param-value pairs.
%
%  FileName = savePlot(Gx)
%
%  FileName = savePlot(Gx, Param, Value, ...)
%
%  INPUT
%    Gx: figure handle to save
%    
%    Parameter    Values(*)   Description
%    ------------ ----------  ---------------------------------------------
%    SaveAs       * ''        Ask user to choose file name.
%                   [str]     File name to save as.
%    SavePre      * ''        Add no prefix.
%                   [str]     Add this string BEFORE the file extension.
%                             See notes below for more details.
%    Format       * 'png'     PNG portable network graphics
%                   'tif'     TIFF image, compressed
%                   'jpg'     JPEG format
%                   'fig'     FIG matlab format
%                   'emf'     EMF enhanced metafile (vector graphics)
%                   'pdf'     PDF document
%    DPI          * 300       Dots per square inch
%                   int       Integer number for DPI. Must be squarable.
%
%  OUTPUT
%    SaveName: the actually file name used to save the image. Can be
%      different from SaveAs if using the SavePre option.
%
%  NOTE
%    SavePre is used when generating many images with similar file names,
%    except a small incremental change. For example, Tree1.png, Tree2.png,
%    Tree3.png, etc. To use save pre, specify the SaveAs name AND SavePre
%    number. For example, savePlot(gcf, 'SaveAs', 'tree.png', SavePre, 1).
%
%  EXAMPLE
%    Gx = loadDemoPlot;
%    savePlot(Gx, 'SaveAs', '', 'DPI', 600);

function varargout = savePlot(varargin)
[Gx, varargin] = getOpenGCF(varargin{:});
ValidFmt = {'tif', 'jpg', 'png', 'fig', 'emf', 'pdf', 'eps', 'svg', 'ps'};
P = inputParser;
addParameter(P, 'SaveAs', '', @ischar);
addParameter(P, 'SavePre', '', @ischar);
addParameter(P, 'Format', 'png', @(x) any(validatestring(lower(x), ValidFmt)));
addParameter(P, 'DPI', 300, @isnumeric);
[Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInput(P, varargin{:});
if ReturnThis
   varargout = {Ps, Pu, ExpPs, ExpPu};
   return
end
SaveAs = Ps.SaveAs;
SavePre = Ps.SavePre;
Format = Ps.Format;
DPI = Ps.DPI;

%Determine appropriate file path, name, and extension
if isempty(SaveAs) %Need to ask users to select file name
    [SaveFile, SavePath] = uiputfile(sprintf('*.%s;', ValidFmt{:}));
    if isnumeric(SaveFile); return; end
    [SavePath, SaveFile, SaveExt] = parseFileName(fullfile(SavePath, SaveFile));
else
    [SavePath, SaveFile, SaveExt] = parseFileName(SaveAs);
    if isempty(SaveExt) || ~ismember(lower(SaveExt), cellfun(@(x) sprintf('.%s', x), ValidFmt, 'un', 0))
        SaveExt = ['.' Format]; %Use the format option specified in the input. Format cannot be empty.
        SaveFile = [SaveFile SaveExt];
    end
end

%Remove the double dot '..jpg' problems
CurLen = length(SaveFile);
SaveFile = strrep(SaveFile, '..', '.');
while CurLen ~= length(SaveFile)
    CurLen = length(SaveFile);
    SaveFile = strrep(SaveFile, '..', '.'); 
end

%In case you have to add a prefix, adjust the file name
if ~isempty(SavePre)
    DotLoc = find(SaveFile == '.');
    SaveFile = [SaveFile(1:DotLoc(end)-1) SavePre SaveFile(DotLoc(end):end)];
end

%Save the file
FullSaveName = fullfile(SavePath, SaveFile);
switch SaveExt  
    case '.tif'
        print(Gx, FullSaveName, '-dtiff', ['-r' num2str(DPI)], '-painters'); %NOTE: painters renderer prevents plot cutoff bugs in matlab.
    case '.jpg'
        print(Gx, FullSaveName, '-djpeg', ['-r' num2str(DPI)], '-painters');
    case '.png'
        print(Gx, FullSaveName, '-dpng', ['-r' num2str(DPI)], '-painters');
    case '.pdf'
        print(Gx, FullSaveName, '-dpdf');
    case '.eps'
        print(Gx, FullSaveName, '-deps');
    case '.svg'
        print(Gx, FullSaveName, '-dsvg');
    case '.ps'
        print(Gx, FullSaveName, '-dps');
    case '.emf'
        saveas(Gx, FullSaveName);
    case '.fig'
        saveas(Gx, FullSaveName);
end

varargout{1} = FullSaveName;