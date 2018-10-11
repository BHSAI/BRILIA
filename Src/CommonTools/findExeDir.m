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
%     if isunix
%         [Status, PsMsg] = system('ps -C BRILIA');
%         assert(Status == 0, '%s: Error with using "pc -C BRILIA."', mfilename);
%         ProcID = regexp(PsMsg, '\d+', 'match', 'once');
%         [Status, LsMsg] = system(sprintf('ls -l /proc/%s/exe', ProcID));
%         assert(Status == 0, '%s: Error with using "ls -l /proc/ID/exe."', mfilename);
%         ExeDir = char(regexpi(LsMsg, '->\s*(.*)', 'tokens', 'once'));
%     elseif ispc
%         [Status, SetMsg] = system('set PATH');
%         assert(Status == 0, '%s: Error with using "set PATH".', mfilename);
%         ExeDir = char(regexpi(SetMsg, 'Path=(.*?);', 'tokens', 'once'));
%     elseif ismac
%         error('%s: Unsupported operating system - mac.', mfilename);
%     else
%         error('%s: Unknown or unsupported operating system.', mfilename);
%     end
else
    ExeDir = findRoot();
end

if nargin == 1 && strcmpi(varargin{1}, 'print')
    fprintf('%s\n', ExeDir);
end