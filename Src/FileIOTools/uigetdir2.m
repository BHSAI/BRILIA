%uigetdir2 will open multiple directories using JAVA's file chooser.
%
%  PathNames = uigetdir2(StartPath, Title)
%
%  INPUT
%    StartPath: the initial path to start the dialog. Defaults to pwd.
%    DialogTitle: the title for the dialog. Defaults to 'Open directories'.
%
%  OUTPUT
%    PathNames: path names to multiple directories, stored in cell array.
%
%Modified from:
%https://www.mathworks.com/matlabcentral/fileexchange/32555-uigetfile-n-dir---select-multiple-files-and-directories
function PathNames = uigetdir2(StartPath, DialogTitle)
import javax.swing.JFileChooser;
if nargin == 0 || isempty(StartPath) || StartPath == 0 % Allow a null argument.
    StartPath = pwd;
end

jchooser = javaObjectEDT('javax.swing.JFileChooser', StartPath);
jchooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
if nargin > 1
    jchooser.setDialogTitle(DialogTitle);
end
jchooser.setMultiSelectionEnabled(true);
status = jchooser.showOpenDialog([]);

if status == JFileChooser.APPROVE_OPTION
    jFile = jchooser.getSelectedFiles();
	PathNames{size(jFile, 1)}=[];
    for i = 1:size(jFile, 1)
        PathNames{i} = char(jFile(i).getAbsolutePath);
    end
elseif status == JFileChooser.CANCEL_OPTION
    PathNames = [];
else
    error('Error occured while picking file.');
end
