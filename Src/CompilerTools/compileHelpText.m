%compileHelpText will search through all m files and directories, pull out
%the help text, and save it to a \HelpText\mfilename.txt. This is used with
%showHelp(mfilename) to show the help file in a command line environment.

function compileHelpText(varargin)
P = inputParser;
addParameter(P, 'Dir', '', @ischar);
addParameter(P, 'CheckSub', 'y', @(x) ischar(x) && ismember(lower(x), {'y', 'n'}));
addParameter(P, 'SaveTo', '', @ischar);
addParameter(P, 'Overwrite', 'ask', @(x) ischar(x) && ismember(lower(x), {'y', 'n', 'ask'}));
parse(P, varargin{:});
P = P.Results;

%Determine if to overwrite existing help files
switch lower(P.Overwrite)
    case 'y'
        OverwriteAll = 1;
    case 'n'
        OverwriteAll = 0;
    otherwise %Ask user only if there's a overwrite
        OverwriteAll = -1;
end

%Format the directory name, to ensure it ends with slash
if isempty(P.Dir)
    P.Dir = cd;
end
if ~strcmp(P.Dir(end), filesep)
    P.Dir = [cd filesep];
end

%Determine the save directory
if isempty(P.SaveTo)
    P.SaveTo = [P.Dir 'HelpText' filesep];
else
    [FilePath, FileName, FileExt] = parseFileName(P.SaveTo, 'ignorefilecheck');
    if isempty(FileExt) 
        if ~isempty(FileName)
            FilePath = cat(2, FilePath, FileName, filesep);
        end
    else
        error('%s: SaveTo should be a folder, not a file name', mfilename);
    end
    P.SaveTo = FilePath;
end
[Success, Msg] = mkdir(P.SaveTo);
assert(Success > 0, '%s: Could not create HelpText folder at "%s".\n  %s', mfilename, TargetDir, Msg)

%Collect subfolder folder into P.Dir, which is now a cell array of folders
if strcmpi(P.CheckSub(1), 'y')
    SubDir = genpath(P.Dir);
    P.Dir = regexp(SubDir, pathsep, 'split')';
    if isempty(P.Dir{end})
        P.Dir(end) = [];
    end
    %Add the '\'
    for j = 1:length(P.Dir)
        if ~strcmp(P.Dir{j}(end), filesep)
            P.Dir{j} = [P.Dir{j} filesep];
        end
    end
else
    P.Dir = {P.Dir};
end

%For every .m file, take the comments before the function and save it into
%a HelpText folder with the same name as the .m file.
for j = 1:length(P.Dir)
    MFile = dir([P.Dir{j} '*.m']);
    for k = 1:length(MFile)
        MFileName = MFile(k).name;
        HelpTextName = [MFileName(1:end-1) 'txt'];
        FileExist = exist([P.SaveTo HelpTextName], 'file');
        if FileExist && OverwriteAll == -1
            AskText = sprintf('%s exists. Do you want to overwrite all? y or n. ', HelpTextName);
            while OverwriteAll == -1
                AskOverwriteAll = input(AskText, 's');
                if isempty(AskOverwriteAll); continue; end
                if strcmpi(AskOverwriteAll(1), 'y')
                    OverwriteAll = 1;
                    break;
                elseif strcmpi(AskOverwriteAll(1), 'n')
                    OverwriteAll = 0;
                    break;
                end
            end
        elseif FileExist && OverwriteAll == 0
            continue;
        end
        
        %Go and write new help text file
        FID1 = fopen([P.Dir{j} MFile(k).name], 'r');
        if FID1 <= 0
            warning('%s: Could not read [ %s ]', mfilename, [P.Dir{j} MFileName]);
        end
        FID2 = fopen([P.SaveTo HelpTextName], 'w');
        if FID2 <= 0
            warning('%s: Could not create [ %s ]', mfilename, [P.Dir{j} HelpTextName]);
            fclose(FID1);
        end
        TextLine = fgetl(FID1);
        while ischar(TextLine)
            StartPos = regexpi(TextLine, '\S', 'once');
            if ~isempty(StartPos)
                if TextLine(StartPos) == '%'
                    fprintf(FID2, '%s\r\n', TextLine(2:end));
                elseif ~isempty(TextLine)
                    EndPos = StartPos + length('function') - 1;
                    if EndPos <= length(TextLine) && strcmpi(TextLine(StartPos:EndPos), 'function')
                        break
                    end
                end
            end
            TextLine = fgetl(FID1);
        end
        fclose(FID1);
        fclose(FID2);
    end
end

