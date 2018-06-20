%findRoot will find the root directory for BRILIA. Returns the
%directory name, ending with the forward or reverse slash, for BRILIA
%starting root. The root path is found relative to the this function,
%by looking for where "Src" or "BRILIAvX.Y.Z" folder is.
%
%  RootDir = findRoot
%
function RootDir = findRoot
FilePath = fileparts(mfilename('fullpath'));
PathParts = strsplit(FilePath, filesep);
MainIdx = find(cellfun(@(x)~isempty(regexpi(x, 'BRILIA.\d+.\d+.\d+')), PathParts), 1, 'last');
SrcIdx = find(strcmp(PathParts, 'Src'), 1, 'last');
if ~isempty(MainIdx)
    RootDir = fullfile(PathParts{1:MainIdx(end)}, filesep); 
elseif ~isempty(SrcIdx)
    RootDir = fullfile(PathParts{1:SrcIdx(end)-1}, filesep);
else
    error('%s: Could not find BRILIA root dir relative to "%s".', mfilename, FilePath);
end