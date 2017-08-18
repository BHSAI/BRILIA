%findRoot will find the root directory for BRILIA. Returns the
%directory name, ending with the forward or reverse slash, for BRILIA
%starting root. The root path is found relative to the this function,
%by looking for where "Src" or "BRILIAvX.Y.Z" folder is.
%
%  RootDir = findRoot()
function RootDir = findRoot
%Get the file name and segment by slash locations
[FilePath, ~, ~] = parseFileName(mfilename('fullpath'));
SlashType = FilePath(regexp(FilePath, '\\|\/', 'once'));
SlashLoc = regexp(FilePath, '\\|\/');
FilePathParts = regexp(FilePath, SlashType, 'split');

%Identify the locations of either the main folder or Src folder (in case
%the user renamed the folder).
MainPart = 0;
SrcPart = 0;
for j = length(FilePathParts):-1:1
    if ~isempty(regexpi(FilePathParts{j}, 'Src', 'once'))
        SrcPart = j;
    end
    if ~isempty(regexpi(FilePathParts{j}, 'BRILIA.\d+.\d+.\d+', 'once'))
        MainPart = j;
        break
    end
end
if MainPart == 0 && SrcPart > 0
    MainPart = SrcPart - 1;
end
if MainPart <= 0
    error('%s: Could not find the root BRILIA folder relative to %s.\nMake sure you do not change the structure of the BRILIA folder,\nas it looks for folder names "BRILIAvX.Y.Z" or "Src".', mfilename, FilePath);
end
RootDir = FilePath(1:SlashLoc(MainPart));
