%prepSaveTarget will prepare the folder and/or file name to use for
%saving items. This was created to simplify making decisions for
%where to save a file given a certain combination of inputs. This
%function will also attempt to create the directory to save the target
%files.
%
%  [FullSaveName, ValidTarget] = prepSaveTarget(Param, Value, ...)
%
%  INPUT
%    Param-Value pairs are as follows:
%    Param        Value    Description
%    -----------  -------  ------------------------------------------------
%    SaveAs         ''     Full or just file name to use w or w/o extension
%    SavePrefix     ''     Any string you want to add between the save file
%                            name and the extension (EX: If SavePrefix =
%                            '.grp1', then save.png => save.grp1.png
%    SaveExt        ''     The file extension to use. See note.
%    SaveDir        ''     Save directory to use for saving files. Will use
%                            current dir if SaveAs has not file path and
%                            SaveDir is empty.
%    SaveSubDir     ''     Sub directory to add on top of any dir in SaveAs
%                            or DefaultDir. Used to simplify adding a sub
%                            dir under an original file dir.
%    MakeSaveDir    'n'    Will only return the FullSaveName.
%                   'y'    Will attempt to create the save dir. 
%    DefaultSaveAs  ''     Will use this if SaveAs is empty. Hint: Let the
%                            function that is going to actually write the
%                            file specify as SaveSubDir and DefaultSaveAs.
%
%  OUTPUT
%    FullSaveName: The full file path and name to use for saving the file
%
%  NOTE
%    SaveExt does not have to be specified if SaveAs has the extension, BUT
%    if specified, will override what is in SaveAs extension. This is to
%    allow users to use a source file name for the SaveAs value, and to
%    simply change the extension to an image. 
%
%    The default order of preference for saving are as follows:
%    SavePath = SaveAs folder > SaveDir folder > or DefaultDir (cd)
%    *Having a SaveSubDir will ALWAYS created a subdir
%    SaveName = SaveAs file name > error
%    SaveExt = This overrides SaveAs' file ext > SaveAs file ext > error
%    
%  EXAMPLE
%    FileName = 'c:\local\file.csv'
%    FullSaveName = prepSaveTarget('SaveAs', FileName, 'SaveExt', ...
%                                               '.png', 'SavePrefix', '.2')
%    FullSaveName = 
%          'c:\local\file.2.png'
%
%    SaveName = 'file.csv'
%    SavePre = '.3'
%    SavePath = 'c:\local'
%    SaveSubDir = '\Plots\New'
%    FullSaveName = prepSaveTarget('SaveAs', SaveName, 'SaveSubDir', ...
%                    SaveSubDir, SavePrefix', SavePre, 'SaveDir', SavePath)
%    FullSaveName = 
%          'c:\local\Plots\file.3.csv'

function varargout = prepSaveTarget(varargin)
P = inputParser;
addParameter(P, 'SaveAs', '', @(x) ischar(x) || isempty(x));
addParameter(P, 'SavePrefix', '', @(x) ischar(x) || isempty(x));
addParameter(P, 'SaveExt', '', @(x) ischar(x) || isempty(x));
addParameter(P, 'SaveDir', '', @(x) ischar(x) || isempty(x));
addParameter(P, 'SaveSubDir', '', @(x) ischar(x) || isempty(x));
addParameter(P, 'DefaultSaveAs', '', @(x) ischar(x) || isempty(x));
addParameter(P, 'MakeSaveDir', 'n', @(x) ischar(x) && ismember(lower(x(1)), {'y', 'n'}));
[Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInput(P, varargin{:});
if ReturnThis
   varargout = {Ps, Pu, ExpPs, ExpPu};
   return;
end
SaveAs = Ps.SaveAs;
SavePrefix = Ps.SavePrefix;
SaveExt = Ps.SaveExt;
SaveDir = Ps.SaveDir;
MakeSaveDir = Ps.MakeSaveDir;
SaveSubDir = Ps.SaveSubDir;
DefaultSaveAs = Ps.DefaultSaveAs;

%Begin with SaveAs
if isempty(SaveAs)
    if isempty(DefaultSaveAs)
        error('%s: Must specify a SaveAs file name.');
    else %Take from the default
        SaveAs = DefaultSaveAs;
    end
end
[FilePath, FileName, FileExt] = parseFileName(SaveAs, 'ignorefilecheck'); %Note that FilePath will return curdir if SaveAs does not have the filepaths in it
if isempty(FileName) 
    error('%s: The SaveAs name [ %s ] is not valid', mfilename, SaveAs);
end
SlashType = FilePath(regexpi(FilePath, '\\|\/', 'once'));

%Override the file path IF SaveAs does not have a file path slash
if ~isempty(SaveDir) && isempty(regexpi(SaveAs, '\\|\/', 'once'))
    FilePath = SaveDir;
    if FilePath(end) ~= SlashType
        FilePath(end+1) = SlashType;
    end
end

%Remove the FileExt from FileName
ExtLoc = regexpi(FileName, ['\' FileExt]);
if ~isempty(ExtLoc)
    FileNamePre = FileName(1:ExtLoc(end)-1);
else
    FileNamePre = FileName;
end

%Check the file extension override
if ~isempty(SaveExt)
    FileExt = SaveExt;
else
    if isempty(FileExt)
        error('%s: Cannot decide the file extension with SaveAs [ %s ] and SaveExt [ %s ] inputs.', mfilename, SaveAs, SaveExt);
    end
end

%Add the subdir, but make sure SaveSubDir is valid
if ~isempty(SaveSubDir)
    %Can't have a double slash
    if ~isempty(regexpi(SaveSubDir, '\\{2,Inf}|\/{2,Inf}', 'once'))
        error('%s: SaveSubDir [ %s ]cannot have double slashes.', mfilename, SaveSubDir);
    end
    %Convert wrong-direction slashes and remove 1st and last slashes
    SlashLoc = regexpi(SaveSubDir, '\\|\/');
    SaveSubDir(SlashLoc) = SlashType;
    FirstLastLoc = intersect(SlashLoc, [1 length(SaveSubDir)]);
    SaveSubDir(FirstLastLoc) = [];
    %Remove invalid char, which are non-letter, non-underscore, non-slash
    InvalidChar = regexpi(SaveSubDir, '[^\w\_\\\/]');
    if ~isempty(InvalidChar)
        warning('%s: Removing characters from SaveSubDir name [ %s ]', mfilename, SaveSubDir(InvalidChar));
        SaveSubDir(InvalidChar) = [];
        warning('%s: SaveSubDir name is [ %s ]', mfilename, SaveSubDir);
    end
    if ~isempty(SaveSubDir)
        FilePath = [FilePath SaveSubDir SlashType];
    end
end

if strcmpi(MakeSaveDir, 'y')
    [Success, Msg] = mkdir(FilePath);
    if Success == 0
        error('%s: Could not create save directory [ %s ].\n  %s', mfilename, FilePath, Msg);
    end
end

varargout{1} = [FilePath FileNamePre SavePrefix FileExt];
