%findRoot will find the root directory for BRILIA. Returns the
%directory name, ending with the forward or reverse slash, for BRILIA
%starting root. The root path is found relative to the this function,
%by looking for where "Src" or "BRILIAvX.Y.Z" folder is.
%
%  RootDir = findRoot
function RootDir = findRoot
FilePath = fileparts(mfilename('fullpath'));
PathParts = regexp(FilePath, filesep, 'split');
SrcLoc = find(cellfun(@(x)strcmp(x, 'Src'), PathParts));
MainLoc = cellfun(@(x)~isempty(regexpi(x, 'BRILIA.\d+.\d+.\d+')), PathParts);
if ~isempty(SrcLoc)
    RootDir = fullfile(PathParts{1:SrcLoc(end)-1}, filesep);
elseif ~isempty(MainLoc)
    RootDir = fullfile(PathParts{1:MainLoc(end)}, filesep); 
else
    error('%s: Could not find BRILIA root dir relative to "%s".', mfilename, FilePath);
end
