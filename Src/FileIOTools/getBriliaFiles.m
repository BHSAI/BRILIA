%getBriliaFiles will search through either a directory or file for the
%final BRILIA annotation file, which should end like *.BRILIAvN*.csv. 
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

function FileNames = getBriliaFiles(FileNames, MultiOn, CmdOn)
if nargin < 1
    FileNames = '';
end
if nargin < 2
    MultiOn = 1;
end
if nargin < 3
    CmdOn = 1;
end

if isempty(dir2(FileNames))
    FileNames = uigetdir2(FileNames, MultiOn, CmdOn, 'Select the BRILIA output files or folders.');
elseif ischar(FileNames)
    FileNames = dir2(FileNames); %{FileNames};
end

for f = 1:length(FileNames)
    if isdir(FileNames{f})
        %Try level 1 subfolder
        TmpFiles = dir2(fullfile(FileNames{f}, '*BRILIAv*.csv'), 'file');
        %Try level 2 subfolder (no more than this)
        if isempty(TmpFiles)
            TmpFiles = dir2(fullfile(FileNames{f}, '*/*BRILIAv*.csv'), 'file');
        end
        GoodLoc = ~endsWith(TmpFiles, {'Raw.csv', 'Err.csv'}, 'ignorecase', true);
    else 
        TmpFiles = dir2(FileNames{f}); %Wrap it up for simply vertcat at end. Use dir2 to check for file existence, which returns empty for non-existing file.
        GoodLoc = ~cellfun('isempty', regexp(TmpFiles, 'BRILIAv\w+.csv'));
    end
    if sum(GoodLoc) > 1 %Perhaps there is a modification
        GoodIdx = find(GoodLoc);
        fprintf('%s: Found multiple BRILIAv*.csv files. Picking the longest-name file.\n', mfilename);
        FileLen = cellfun('length', TmpFiles(GoodIdx));
        [~, KeepIdx] = max(FileLen);
        FileNames{f} = TmpFiles(GoodIdx(KeepIdx(1)));
    else
        FileNames{f} = TmpFiles(GoodLoc);
    end
end
FileNames = vertcat(FileNames{:});
