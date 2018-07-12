%savePlot will save the plot specified by the gcf, save name, and save
%param-value pairs.
%
%  FileName = savePlot(Gx)
%
%  FileName = savePlot(Gx, Param, Value, ...)
%
%  INPUT
%    Gx: figure handle to save
%    P: strucutre of default input values
%    'getinput': returns only the default inputs of savePlot as a structure
%    'getvarargin': returns param-value pairs used by savePlot
%    varargin: param-value pairs that are used by savePlot and others. Use
%      'getvarargin' to remove extra param-value pairs prior to running
%      savePlot(varargin).
%    
%    Param-Value pairs are defined as:   
%      Parmaeter Name  Value        Description
%      --------------- ---------    ---------------------------------------
%      Save            'n', 'y'     To save or not to save. Default 'y'.
%      SaveAs          String       The name to save the file. Will add
%                                     appropriate extensions.
%      SavePre         String       Any additional name to add BEFORE the
%                                     extension. See notes.
%      Format          'tif' 'png'  Image format type to save  figure as.
%                      'jpg' 'fig'    Default is 'tif'.
%      DPI             Integer      Dots per square inch. Default is 400,
%                                     which must be a square number.
%
%  OUTPUT
%    P: structure of default inputs of this function
%    SaveName: the actually file name used to save the image. Can be
%      different from SaveAs due to modifiers and corrections.
%    SaveVarargin: varargin used only by savePlot.
%    varargin: remaining varargin parameters not used by savePlot.
%
%  NOTE
%    SavePre is used when generating many images with the same file name,
%    but incremental changes. Ex: 'SaveAs' Tree.png, but want Tree1.png,
%    Tree2.png, Tree3.png... In this case, 'SavePre' is used as '1', '2',
%    '3'... per each savePlot summon. It's used mainly to simplify having
%    to parse and remake file names.

function varargout = savePlot(varargin)
[Gx, varargin] = getOpenGCF(varargin{:});
P = inputParser;
addParameter(P, 'Save', 'y', @(x) any(validatestring(lower(x), {'y', 'n'})));
addParameter(P, 'SaveAs', '', @ischar);
addParameter(P, 'SavePre', '', @ischar);
addParameter(P, 'Format', 'png', @(x) any(validatestring(lower(x), {'tif', 'jpg', 'png', 'fig'})));
addParameter(P, 'DPI', 300, @isnumeric);
[Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInput(P, varargin{:});
if ReturnThis
   varargout = {Ps, Pu, ExpPs, ExpPu};
   return;
end
Save = Ps.Save;
SaveAs = Ps.SaveAs;
SavePre = Ps.SavePre;
Format = Ps.Format;
DPI = Ps.DPI;

%Save figures if needed
if strcmpi(Save, 'y')
    %Determine appropriate file path, name, and extension
    if isempty(SaveAs) %Need to ask users to select file name
        [SaveFile, SavePath] = uiputfile('*.tif;*.png;*.jpg;*.fig');
        if isnumeric(SaveFile)
            return;
        end
        [SavePath, SaveFile, SaveExt] = parseFileName([SavePath SaveFile], 'ignorefilecheck');
    else
        [SavePath, SaveFile, SaveExt] = parseFileName(SaveAs, 'ignorefilecheck');
        if isempty(SaveExt) || ~ismember(lower(SaveExt), {'.tif', '.png', '.jpg', '.fig'})
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
    
    %Assemble the full save name, and save depending on file ext.
    %NOTE: Use the painters renderer to prevent plot cutoff bugs in matlab.
    FullSaveName = fullfile(SavePath, SaveFile);
    switch SaveExt
        case '.tif'
            print(Gx, FullSaveName, '-dtiff', ['-r' num2str(DPI)], '-painters');
        case '.jpg'
            print(Gx, FullSaveName, '-djpeg', ['-r' num2str(DPI)], '-painters');
        case '.png'
            print(Gx, FullSaveName, '-dpng', ['-r' num2str(DPI)], '-painters');
        case '.fig'
            saveas(Gx, FullSaveName);
    end
    
    if nargout >= 1
        varargout{1} = FullSaveName;
    end
end
