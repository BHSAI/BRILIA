%dir2 is a form of dir that returns the full file name and enables
%filtering out the '.' and '..' results from normal dir. 
%
%  FullName = dir2(FileStr)
%
%  FullName = dir2(FileStr, Filter)
%
%  [FullName, NameOnly, PathOnly] = dir2(FileStr, ...)
%
%  INPUT
%    FileStr: file name or string to place into dir2
%    Filter ['file', 'dir']: select only files or dir
%    
%  OUTPUT
%    FullName: cell of full directory or file names
%    NameOnly: cell of only the name of the directories or files
%    PathOnly: cell of only the path before the directories or files
%
function [FullName, NameOnly, PathOnly] = dir2(FileStr, Filter)
if nargin == 1
    Filter = 'none';
elseif nargin == 2 && ~ismember(lower(Filter), {'file', 'dir'})
    error('%s: Filter must be ''file'', or ''dir''.', mfilename);
end

DirOut = dir(FileStr);
DotLoc = ismember({DirOut.name}, {'.', '..'});
DirOut = DirOut(~DotLoc);

switch lower(Filter)
    case 'file'
        DirOut = DirOut(~[DirOut.isdir]);
    case 'dir'
        DirOut = DirOut([DirOut.isdir]);
end

FullName = fullfile({DirOut.folder}, {DirOut.name})';
if nargout >= 2
    NameOnly = {DirOut.name}';
    if nargout >= 3
        PathOnly = {DirOut.folder}';
    end
end