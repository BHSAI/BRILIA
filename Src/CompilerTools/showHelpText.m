%showHelpText will display a m file's help text on the command line, which
%should be stored in a 'HelpText' folder. This is a workaround solution for
%matlab's inability to show help text in compiled codes.
%
%  showHelpText('Name', Name, 'HelpTextDir', HelpTextDir)
%
%  INPUT
%    Name: the name of the function or m file (without the .m)
%    HelpDir: the directory storing all help text files created from
%      compileHelpText
%
%  See also compileHelpText

function showHelpText(Name, varargin)
CurDir = cd;
SlashType = CurDir(regexpi(CurDir, '\\|\/', 'once'));

%If HelpDir is given, see if it exists
HelpDir = '';
if ~isempty(varargin) && ischar(varargin{1})
    HelpDir = varargin{1};
    if HelpDir(end) ~= SlashType
        HelpDir = cat(2, HelpDir, SlashType);
    end
    if ~exist(HelpDir, 'dir');
        HelpDir = '';
    end
end

%Determine if the file exists (without path)
DotLoc = find(Name == '.');
if ~isempty(DotLoc) %See if it ends in .m
    if ismember(lower(Name(DotLoc(end):end)), {'.m', '.txt'})
        Name = Name(1:DotLoc(end) - 1);
    end
end
HelpTextName = [Name '.txt'];

if ~exist([HelpDir HelpTextName], 'file') %Need to search for a valid HelpDir
    %Look for the HelpDir
    SubDir = genpath(CurDir);
    SubDir = regexp(SubDir, ';', 'split')';
    if isempty(SubDir{end})
        SubDir(end) = [];
    end
    
    HelpDirLoc = zeros(length(SubDir), 1);
    for j = 1:length(SubDir)
        FolderParse = regexp(SubDir{j}, '\\|\/', 'split');
        for k = length(FolderParse):-1:1
            if ~isempty(regexpi(FolderParse{k}, 'HelpText', 'once'))
                HelpDirLoc(j, 1) = k;
                break;
            end
        end
    end
    
    PotentialLoc = find(HelpDirLoc > 1);
    if isempty(PotentialLoc)
        return;
    end
    
    SelectLoc = PotentialLoc(HelpDirLoc(PotentialLoc) == min(HelpDirLoc(PotentialLoc)));
    SelectLoc = SelectLoc(1);
    HelpDir = SubDir{SelectLoc};
    if HelpDir(end) ~= SlashType
        HelpDir = cat(2, HelpDir, SlashType);
    end
end

if ~exist([HelpDir HelpTextName], 'file')
    warning('%s: Could not find help file [ %s ]', mfilename, [HelpDir HelpTextName]);
    return;
end
FID1 = fopen([HelpDir HelpTextName], 'r');
if FID1 <= 0
    warning('%s: Could not read [ %s ]', mfilename, [HelpDir HelpTextName]);
    return
end
TextLine = fgetl(FID1);
while ischar(TextLine)
    fprintf('%s\n', TextLine);
    TextLine = fgetl(FID1);
end
fclose(FID1);
