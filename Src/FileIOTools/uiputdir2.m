%uiputdir2 will ask user where to save the file, either via a UI or command
%line option.
%
%  PathName = uiputdir2(StartPath, DialogTitle)
%
%  PathName = uiputdir2(..., 'cmd')
%
%  INPUT
%    StartPath: the initial path to start the dialog. Defaults to pwd.
%    DialogTitle: the title for the dialog box. 
%    'cmd': User must use command line to select the folder to use.
%
%  OUTPUT
%    PathName: full path to the output directory
%
%  NOTE
%    If JAVA is disabled, default option becomes 'cmd'.
%
%  EXAMPLE
%    DirName = uiputdir2(pwd, 'Select Dir');
%
%    DirName = uiputdir2(pwd, 'Select Dir', 'cmd');
%
function PathName = uiputdir2(varargin)
persistent StartPath   %The most recent folder is persistent
if isempty(StartPath)
    StartPath = pwd;
end

%Check for cmd option
CmdLoc = ismember(lower(varargin), 'cmd');
CmdOn = any(CmdLoc);
varargin = varargin(~CmdLoc);

%Establish first search folder 
if numel(varargin) >= 1 && isdir(varargin{1})
    StartPath = varargin{1};
end

%Set the message prompt
if numel(varargin) >= 2
    DialogTitle = varargin{2};
else
    DialogTitle = 'Choose directory to save to';
end

%Find valid folder(s)
PathName = {};
if usejava('jvm') && ~CmdOn   %For using the GUI
    import javax.swing.JFileChooser;
    import javax.swing.filechooser.FileNameExtensionFilter;

    JFileChooser = javaObjectEDT('javax.swing.JFileChooser', StartPath);
    JFileChooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
    JFileChooser.setMultiSelectionEnabled(false);
    JFileChooser.setDialogTitle(DialogTitle);

    Status = JFileChooser.showOpenDialog([]);
    if Status == JFileChooser.APPROVE_OPTION
        JFile = JFileChooser.getSelectedFile();
        PathName = char(JFile.getAbsolutePath);
        if ~isdir(PathName)
            PathName = fileparts(PathName);
        end
        StartPath = PathName;
    elseif Status == JFileChooser.CANCEL_OPTION
        fprintf('Cancelled folder selection.\n');
    end
    
else %Ask user to type the directory path
    while isempty(PathName)
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

        else %User has specified a string for folders
            if ispc && ~contains(PathInput, ':') && PathInput(1) ~= '.' %User did not specify HDD name
                PathInput = fullfile(['.' filesep], PathInput);         %Append "./" current dir to beginning
            elseif ~startsWith(PathInput, {'.', filesep})
                PathInput = fullfile(['.' filesep], PathInput);         %Append "./" current dir to beginning
            end
            
            PathName = regexprep(PathInput, {'''|"', ['\.\.' filesep], ['\.' filesep]}, strrep({'', [fileparts(StartPath) filesep], [StartPath filesep]}, '\', '\\')); %relative pathing ../ and ./
            if isempty(PathName)
                fprintf('Empty char or invalid folder name was given.\n');
            end
        end
    end
end