%uigetdir2 will open multiple directories using JAVA's file chooser.
%
%  PathNames = uigetdir2(StartPath, Title)
%
%  PathNames = uigetdir2(StartPath, Title, FilterExt)
%
%  INPUT
%    StartPath: the initial path to start the dialog. Defaults to pwd.
%    DialogTitle: the title for the dialog. Defaults to 'Open directories'.
%    FilterExt: a cell array of filter extensions to use. Ex: {'csv','txt'}
%
%  OUTPUT
%    PathNames: path names to multiple directories, stored in cell array.
%
%Modified from:
%https://www.mathworks.com/matlabcentral/fileexchange/32555-uigetfile-n-dir---select-multiple-files-and-directories
function PathNames = uigetdir2(StartPath, DialogTitle, varargin)
import javax.swing.JFileChooser;
import javax.swing.filechooser.FileNameExtensionFilter;
if nargin == 0 || isempty(StartPath) || StartPath == 0 % Allow a null argument.
    StartPath = pwd;
end

%Determine if multi word is used.
MultiLoc = startsWith(varargin, 'multiselect', 'ignorecase', true);
if any(MultiLoc)
    MultiOn = true;
    varargin(MultiLoc) = [];
else
    MultiOn = false;
end

%Instantiate JFileChooser
jchooser = javaObjectEDT('javax.swing.JFileChooser', StartPath);
jchooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
jchooser.setMultiSelectionEnabled(MultiOn);
if nargin >= 2
    jchooser.setDialogTitle(DialogTitle);
end

%Ensure filter are correct
if ~isempty(varargin)
    FilterExt = strrep(varargin, '.', '');
    FilterDesc = sprintf('%s,', FilterExt{:});
    filter = javaObjectEDT('javax.swing.filechooser.FileNameExtensionFilter', FilterDesc(1:end-1), FilterExt);
    jchooser.setFileFilter(filter);
end

%Choose files and store into a cell array of full paths
status = jchooser.showOpenDialog([]);
if status == JFileChooser.APPROVE_OPTION
    if MultiOn
        jFile = jchooser.getSelectedFiles();
    else
        jFile = jchooser.getSelectedFile();
    end
	PathNames{size(jFile, 1)}=[];
    for i = 1:size(jFile, 1)
        PathNames{i} = char(jFile(i).getAbsolutePath);
    end
elseif status == JFileChooser.CANCEL_OPTION
    PathNames = [];
else
    error('Error occured while picking file.');
end