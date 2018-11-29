%showHelp will display a m file's help text on the command line, which
%should be stored in at RootDir\HelpText\. This is a workaround solution
%for matlab's inability to show help text in compiled codes.
%
%  showHelp(Name)
%
%  showHelp(CLIvarargin)
%
%  INPUT
%    Name: the name of the function or m file (without the .m)
%    CLIvarargin: variable input arguments generated from CLI, in which if
%    user types 'help plotTree' or 'info plotTree', or '-h plotTree', then
%    showHelp will attempt to show the help text.
%
function showHelp(varargin)
for j = 1:length(varargin)
    displayHelpText(varargin{j});
end

function displayHelpText(FuncName)
if ~isdeployed
    HelpDir = fullfile(findRoot, 'HelpText');
    HelpFile = fullfile(HelpDir, [FuncName '.txt']);
else
    HelpFile = [FuncName '.txt'];
end
if ~exist(HelpFile, 'file')
    warning('%s: Could not find the help text for "%s" at "%s".', mfilename, FuncName, HelpFile);
    return
end

FID = fopen(HelpFile, 'r');
TXT = textscan(FID, '%s', 'delimiter', '\n');
fclose(FID);
fprintf('\n');
fprintf('%s\n', TXT{1}{:});
fprintf('\n');