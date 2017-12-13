%getBriliaFiles will search through either a directory or file for the
%final BRILIA annotation file, which should end like *.BRILIAvN.csv. 
%
%  INPUT
%    FileNames: char or cell of file names
%
%  OUTPUT
%    FileNames: correct cell of file names leading to the *.BRILIAvN.csv
%      file.

function FileNames = getBriliaFiles(FileNames, SingleOption)
if nargin == 0 || isempty(FileNames)
    if nargin == 2 && startsWith(SingleOption, 'single', 'ignorecase', true)
        FileNames = uigetdir2('', 'Select the BRILIA output files or folders.');
    else
        FileNames = uigetdir2('', 'Select the BRILIA output files or folders.', 'multiselect');
    end
    if numel(FileNames) == 0
        error('%s: No file(s) selected.', mfilename);
    end
elseif ischar(FileNames)
    FileNames = {FileNames};
end
for f = 1:length(FileNames) 
    if exist(FileNames{f}, 'dir')
        Ver = BRILIA('version');
        EndName = ['*BRILIAv' Ver(1) '.csv'];
        DirList = dir(fullfile(FileNames{f}, EndName));
        if isempty(DirList)
            error('%s: Coult not find the file ending with "%s" in dir "%s".', mfilename, EndName, FileNames{f});
        end
        FileNames{f} = fullfile(FileNames{f}, DirList(1).name);
    end
    if ~exist(FileNames{f}, 'file') && ~exist(FileNames{f}, 'dir')
        error('%s: Could not find the file "%s".', mfilename, FileNames{f});
    end
end