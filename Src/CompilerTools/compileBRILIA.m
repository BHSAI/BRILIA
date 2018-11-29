%compileBRILIA will compile BRILIA and generate binary files appropriate
%for the OS (Windows or Linux) being used.

function compileBRILIA(varargin)
RootDir = findRoot;
if ispc
    OS = 'Win';
else
    OS = 'Linux';
end
TargetDir = fullfile(RootDir, 'Bin', OS);
[Success, Msg] = mkdir(TargetDir);
assert(Success > 0, '%s: Could not create directory "%s".\n  %s', mfilename, TargetDir, Msg)

if ~(nargin == 1 && strcmpi(varargin{1}, 'skip'))
    fprintf('%s: %s\n', mfilename, 'Compiling MEX files.');
    compileMexBRILIA;

    fprintf('%s: %s\n', mfilename, 'Making help texts.');
    makeHelpFile;    
end

fprintf('%s: %s\n', mfilename, 'Compiling BRILIA');
CurDir = pwd; %Note: Use CD and then relative folder paths. Otherwise, the executable file retains THIS computer's folder structure, which isn't good. 
cd(findRoot);
mcc('-m', 'BRILIA.m', ...
    '-a', 'Src', ...
    '-a', 'Tables', ...
    '-a', 'HelpText', ...
    '-a', 'Databases', ...
    '-a', 'Examples', ...
    '-d', TargetDir);
cd(CurDir);

fprintf('%s: %s\n', mfilename, 'Cleaning and zipping BRILIA files');
MatlabFiles = {'requiredMCRProducts.txt', 'mccExcludedFiles.log', 'readme.txt'};
DeleteFiles = fullfile(TargetDir, MatlabFiles);
delete(DeleteFiles{:})

zipBRILIA;

fprintf('%s: %s\n', mfilename, ['Finished compiling BRILIA ' BRILIA('version')]);