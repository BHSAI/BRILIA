%newMFile will make a new .m file for matlab codes with a template format.
%Rhis is used to help standardize the comments.
%
%  newMFile(FileName)
%
%  newMFile(FileName, FilePath)
%
%  Input
%    FileName: the name of the m file. Will append ".m" if not given. Can
%      use absolute path as well. FilePath: path of the m file to create.
%    Will use pwd if not given.
%  
function newMFile(FileName, FilePath)
OverWrite = false;
if nargin == 0
    [FileName, FilePath] = uiputfile('*.m', 'Make new M file');
    if isnumeric(FileName); return; end
    OverWrite = true; %uiputfile will ask to overwrite already. 
end

if nargin < 2 || isempty(FilePath) && ~contains(FileName, filesep)
    FilePath = pwd;
end
if ~endsWith(FileName, '.m', 'ignorecase', 1)
    FileName = [FileName '.m'];
end

FullName = fullfile(FilePath, FileName);

if exist(FullName, 'file') && nargin ~= 0
    AskOverwrite = input(sprintf('Do you want to overwrite existing file "%s"? y or n: ', strrep(FullName, '\', '\\')), 's');
    OverWrite = ~(isempty(AskOverwrite) || strcmpi(AskOverwrite(1), 'n'));
    if ~OverWrite
        fprintf('%s: Aborted creating new M file "%s".\n', mfilename, FileName);
        return
    end
end

if OverWrite && exist(FullName, 'file')
    PreviousRecycleState = recycle;
    recycle('on');
    delete(FullName);
    recycle(PreviousRecycleState);
    fprintf('%s: Moved older "%s" to recycle bin.\n', mfilename, FileName);
end

[~, MName] = fileparts(FullName);

Str = {
'%FUNCTION description'
'%'
'%  [O1, O2] = FUNCTION(I1, I2, Parameter, Value, ...)'
'%'
'%  INPUT'
'%    I1:'
'%    I2:'
'%'
'%    Parameter       Value      Details'
'%    --------------- ---------- -------------------------------------------'
'%'
'%  OUTPUT'
'%    O1:'
'%    O2:'
'%'
'function [O1, O2] = FUNCTION(I1, I2, varargin)'};
Str = strrep(Str, 'FUNCTION', MName);

FID = fopen(FullName, 'w');
assert(FID > 0, '%s: Could not create a new m file name "%s".', mfilename, FullName);
fprintf(FID, '%s\n', Str{:});
fclose(FID);