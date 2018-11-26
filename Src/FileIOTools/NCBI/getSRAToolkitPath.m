%getSRAToolkitPath will return the SRA Toolkit path if it is installed in
%the computer. If not found, will ask user to set the SRAT Toolkit path on
%the the local environment variable.
%
%  ToolkitPath = getSRAToolkitPath
%
%  OUTPUT
%    ToolkitPath: The local system path for the SRA Toolkit
%                 Empty, '', if SRA Toolkit was not found in system path
%                 
function ToolkitPath = getSRAToolkitPath
ToolkitPath = '';
SysPath = strsplit(getenv('PATH'), pathsep);
PathIdx = find(contains(SysPath, 'sratoolkit'));
if ~isempty(PathIdx)
    ToolkitPath = sprintf(['%s' pathsep], SysPath{PathIdx});
end
BinPat = sprintf('sratoolkit.*%sbin', ['\' filesep]);
if ~isempty(ToolkitPath) && isempty(regexpi(ToolkitPath, BinPat, 'once'))
    warning('%s: Did not find the sratoolkit''s bin folder.\nTry setSRAToolkitPath to add all toolkit paths.', mfilename);
end
if nargout == 0 && ~isempty(ToolkitPath)
    fprintf('%s\n', ToolkitPath);
end