%openFileDialog will determine if matlab is running headless or not, and
%decide how to ask the user to select a file to open. Note that this does
%not check validity of file name, so that must be done separately after
%using this function.
%
%  [FileName, FilePath] = openFileDialog()
%
%  [FileName, FilePath] = openFileDialog(varargin_uigetfile{:})
%
%  INPUT
%    varargin_uigetfile: inputs for the uigetfile Matlab function
%
%  OUTPUT
%    FileNames: cell of strings of file names (always a cell!)

function FileNames = openFileDialog(varargin)
if usejava('jvm') && ~feature('ShowFigureWindows')
    fprintf('Type file name to open. To cancel, just press empty.\n');
    TriesLeft = 5;
    while 1
        TriesLeft = TriesLeft - 1;
        if TriesLeft < 0
            fprintf('Could not find file "%s".\n Abort.', FileInput);
            FileNames = {};
            return
        end
        FileInput = input('File name: ', 's');
        if exist(FileInput, 'file') || exist(fullfile(pwd, FileInput), 'file')
            break
        else
            fprintf('Could not find file "%s".\n', FileInput);
            continue
        end
    end
    if ~contains(FileInput, filesep)
        FileNames = {fullfile(pwd, FileInput)};
    else
        FileNames = {FileInput};
    end
else
    [FileNames, FilePath] = uigetfile(varargin{:});
    if isnumeric(FileNames)
        FileNames = {};
        return
    end
    if ischar(FileNames); FileNames = {FileNames}; end
    FileNames = cellfun(@(x) fullfile(FilePath, x), FileNames, 'UniformOutput', false)';
end