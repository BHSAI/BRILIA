%delSRAToolkitPath will delete the sratoolkit paths from the PATH
%environmental variable.
%
%  delSRAToolkitPath
%
%
function delSRAToolkitPath
ToolkitPath = getSRAToolkitPath;
if isempty(ToolkitPath); return; end
ToolkitCell = strsplit(ToolkitPath, pathsep);
CurPathCell = strsplit(getenv('PATH'), pathsep);
KeepLoc = ~ismember(CurPathCell, ToolkitCell);
NewPath = sprintf(['%s' pathsep], CurPathCell{KeepLoc});
setenv('PATH', NewPath);