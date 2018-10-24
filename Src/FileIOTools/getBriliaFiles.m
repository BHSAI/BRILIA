%getBriliaFiles will search through either a directory or file for the
%final BRILIA annotation file, which should end like *.BRILIAvN.csv. 
%
%  FileNames = getBriliaFiles(FileNames)
%
%  FileNames = getBriliaFiles(FileNames, Option)
%
%  INPUT
%    FileNames: cell array of file or folder names. If empty, ask user to
%      select.
%    Option: determine how to ask user for files (from uigetdir2.m)
%      'multiselect' : User can select multiple folders
%      'single': User can select only one folder
%      'cmd' : User must use command line to select file(s)
%
%  OUTPUT
%    FileNames: cell array of file names ending with "*.BRILIAvN.csv".
%
%  NOTE
%    If using a custom output file name of BRILIA, user must specify the
%    file names, not folder names.

function FileNames = getBriliaFiles(FileNames, Option)
if nargin < 2
    Option = 'multiselect';
elseif ~ismember(lower(Option), {'multiselect', 'cmd', 'single'})
    error('%s: Option must be ''multiselect'', ''cmd'', or ''single''.', mfilename);
end   

if nargin == 0 || isempty(FileNames)
    FileNames = uigetdir2('', 'Select the BRILIA output files or folders.', Option);
    if isempty(FileNames)
        warning('%s: No file(s) selected.', mfilename);
        return
    end
elseif ischar(FileNames)
    FileNames = {FileNames};
end
    
for f = 1:length(FileNames) 
    if isdir(FileNames{f})
        FileNames{f} = dir2(fullfile(FileNames{f}, '*BRILIAv*.csv'), 'file');
        GoodLoc = ~endsWith(FileNames{f}, {'Raw.csv', 'Err.csv'}, 'ignorecase', true);
        FileNames{f} = FileNames{f}(GoodLoc);
    else 
        FileNames{f} = dir2(FileNames{f}); %Wrap it up for simply vertcat at end
    end
end
FileNames = vertcat(FileNames{:});
