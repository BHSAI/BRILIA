%setSRAToolkitPath will set the path to the SRA toolkit directory, which
%normally has the name like "sratoolkit.2.8.2-1-win64".
%
%  ToolkitPath = setSRAToolkitPath
%
%  ToolkitPath = setSRAToolkitPath(ToolkitPath)
%
%  INPUT
%    ToolkitPath: If empty, will ask use to input the path.
%                 If specified, will attempt to add this path. 
%  
function setSRAToolkitPath(varargin)
if nargin == 0
    PathInput = input('Where is the SRA toolkit folder (ex: c:\Users\usr1\Desktop\sratoolkit.2.8.2-1-win64) \n  SRA Toolkit Path: ', 's');
elseif ischar(varargin{1})
    PathInput = varargin{1};
else
    error('%s: Input must the SRA toolkit directory (ex: ''c:\\Users\\usr\\sratoolkit.2.8.2-1-win64'')', mfilename);
end
assert(isdir(PathInput), '%s: This is not a valid directory - "%s". Recheck spelling', mfilename, PathInput);

%Ensure to add the sratoolkit main folder, and not a subfolder. User could point to subfolder by accident.
PathSegments = strsplit(PathInput, filesep);
SRADirIdx = find(startsWith(PathSegments, 'sratoolkit', 'ignorecase', true), 1, 'last');
assert(~isempty(SRADirIdx), '%s: Could not find the dir with name starting with "sratoolkit" at "%s".', mfilename, PathInput); 

%Add the Toolkit path to the local user environment variable
ToolkitPath = genpath(fullfile(PathSegments{1:SRADirIdx}));
NewPath = strrep([getenv('PATH') pathsep ToolkitPath pathsep], repelem(pathsep, 1, 2), pathsep);
NewPath = unique(strsplit(NewPath, pathsep), 'stable');
NewPath = NewPath(~cellfun('isempty', NewPath));
NewPath = sprintf(['%s' pathsep], NewPath{:});
setenv('PATH', NewPath);