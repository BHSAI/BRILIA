%findExeDir will find the executable directory that hosts this BRILIA. 
%
%  ExeDir = findExeDir
%
%  findExeDir('print')
%
%  INPUT
%    'print': will show the exe-contaning dir on stdout, the command line
%
%  OUTPUT
%    ExeDir: directory full name for the BRILIA executable file. If this is
%      invoked in MATLAB, it will return the root directory of BRILIA
%      source files.
%
%  NOTE
%    For linux, this can return unpredictable results if multiple version
%    of BRILIA are located in separated directories and are also running.

function ExeDir = findExeDir(varargin)
if isdeployed
    ExeDir = fileparts(getExeLocationMEX);
else
    ExeDir = findRoot();
end

if ~ispc
    ExeDir = fullfile(filesep, ExeDir); %Need the initial "/" for some file function
end

if nargin == 1 && strcmpi(varargin{1}, 'print')
    fprintf('%s\n', ExeDir);
end