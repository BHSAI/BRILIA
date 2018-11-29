%uigetdir2 will open multiple directories using JAVA. If no JAVA is
%enabled, will ask user to input the directory path.
%
%  PathNames = uigetdir2(StartPath, DialogTitle)
%
%  PathNames = uigetdir2(..., Option)
%
%  INPUT
%    StartPath: the initial path to start the dialog. Defaults to pwd.
%    DialogTitle: the title for the dialog box. 
%    Option: determine how to ask user for files
%      'multiselect' : User can select multiple folders
%      'single': User can select only one folder
%      'cmd' : User must use command line to select file(s)
%
%  OUTPUT
%    PathNames: cell of full path names to directories.
%
%  NOTE
%    The use of 'cmd' and 'multiselect' are very special and cannot be used
%    as StartPath or DialogTitle! 
%    Do NOT do this: uigetdir2(pwd, 'cmd', 'cmd');
%
%    If JAVA is disabled, default option becomes 'cmd'.
%
%  EXAMPLE
%    DirNames = uigetdir2(pwd, 'Select Dir', 'multiselect');
%
%    DirNames = uigetdir2(pwd, 'Select Dir', 'cmd');
%    Select Dir or 'cancel': ./MyFile
%
function PathNames = uigetdir2(GetDir, MultiOn, CmdOn, DialogTitle)
persistent StartPath   %The most recent folder is persistent
if isempty(StartPath); StartPath = pwd; end

if nargin == 0 
    GetDir = ''; 
else
    GetDir = regexprep(GetDir, '\\|\/', filesep);
end
if nargin < 2; MultiOn = 0; end
if nargin < 3;   CmdOn = 0; end

%Set the message prompt
if nargin < 4
    if MultiOn
        DialogTitle = 'Choose directories to open';
    else
        DialogTitle = 'Choose directory to open';
    end
end

%Find valid folder(s)
if strcmpi(GetDir, ['.' filesep])  %Just return current folder
    PathNames = {StartPath};
elseif strcmpi(GetDir, ['.' filesep '*'])  %Just return subfolders of current folder
    PathNames = dir2(StartPath, 'dir');
elseif isdir(GetDir) && ~contains(GetDir, '*')
    PathNames = {GetDir};
elseif usejava('jvm') && ~CmdOn   %For using the GUI
    [PathNames, StartPath] = getDirUsingJAVA(StartPath, DialogTitle, MultiOn);
else %Ask user to type the directory path
    [PathNames, StartPath] = getDirUsingCMD(StartPath, DialogTitle);
end

function [PathNames, StartPath] = getDirUsingJAVA(StartPath, DialogTitle, MultiOn)
import javax.swing.JFileChooser;
import javax.swing.filechooser.FileNameExtensionFilter;

JFileChooser = javaObjectEDT('javax.swing.JFileChooser', StartPath);
JFileChooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
JFileChooser.setMultiSelectionEnabled(MultiOn);
JFileChooser.setDialogTitle(DialogTitle);

PathNames = {};
Status = JFileChooser.showOpenDialog([]);
if Status == JFileChooser.APPROVE_OPTION
    if MultiOn
        JFile = JFileChooser.getSelectedFiles();
    else
        JFile = JFileChooser.getSelectedFile();
    end
    PathNames = cell(size(JFile, 1), 1);
    for i = 1:numel(JFile)
        PathNames{i} = char(JFile(i).getAbsolutePath);
    end
    StartPath = fileparts(PathNames{1});
elseif Status == JFileChooser.CANCEL_OPTION
    fprintf('Cancelled folder selection.\n');
end
    
function [PathNames, StartPath] = getDirUsingCMD(StartPath, DialogTitle)
fprintf('pwd = "%s"\n', StartPath);
PathNames = {};
while isempty(PathNames)
    PathInput = input(sprintf('%s or ''cancel'': ', DialogTitle), 's');
    PathInput = regexprep(strtrim(PathInput), '\\|/', filesep); %Cleanup input

    if strcmpi(PathInput, 'cancel'); break; end

    if startsWith(PathInput, {'cd', 'pwd', 'ls', 'dir'}, 'ignorecase', true) %User is using folder navigation commands
        InputVar = cleanCommandLineInput(PathInput); %Cleanup input and double quotes
        switch lower(InputVar{1})
            case 'cd'
                if numel(InputVar) > 1
                    if strcmp(InputVar{2}, '..')
                        StartPath = fileparts(StartPath);
                    else
                        if ~startsWith(InputVar{2}, '.')
                            InputVar{2} = fullfile('./', InputVar{2});
                        end
                        TempPath = regexprep(InputVar{2}, {'''|"', ['\.\.' filesep], ['\.' filesep]}, strrep({'', [fileparts(StartPath) filesep], [StartPath filesep]}, '\', '\\'));
                        if isdir(TempPath)
                            StartPath = TempPath;
                        else
                            fprintf('Not a valid directory "%s".\n', InputVar{2});
                        end
                    end
                    if endsWith(StartPath, '..')
                        StartPath = fileparts(StartPath);
                    end
                end
            case 'pwd'
                fprintf('%s\n', StartPath);
            case 'ls'
                ls(StartPath);
            case 'dir'
                dir(StartPath);
        end
        continue

    else
        if ~any(strcmpi(PathInput, {['.' filesep], '.'})) %It is NOT just this folder, but something else
            if ispc && ~contains(PathInput, ':') && PathInput(1) ~= '.' %User did not specify HDD name
                PathInput = fullfile(['.' filesep], PathInput);         %Append "./" current dir to beginning
            elseif ~startsWith(PathInput, {'.', filesep})
                PathInput = fullfile(['.' filesep], PathInput);         %Append "./" current dir to beginning
            end
        end    
        PathInput = regexprep(PathInput, {'''|"', ['\.\.' filesep], ['\.' filesep]}, strrep({'', [fileparts(StartPath) filesep], [StartPath filesep]}, '\', '\\')); %relative pathing ../ and ./
        PathNames = strsplit(PathInput, pathsep)';
        for f = 1:numel(PathNames)
            if endsWith(PathNames{f}, filesep)
                PathNames{f} = PathNames{f}(1:end-1);
            end
            if contains(PathNames{f}, '*')  %Perform wildcard search 
                PathNames{f} = dir2(PathNames{f}, 'dir'); 
            elseif ~isdir(PathNames{f})     %Delete folder that don't exist
                PathNames{f} = {};
            end
        end
        PathNames = vertcat(PathNames{:});
        if ischar(PathNames)
            PathNames = {PathNames};
        elseif isempty(PathNames)
            fprintf('No valid folder selection.\n');
        end
    end
end