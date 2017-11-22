%compileBRILIA will compile BRILIA and generate binary files appropriate
%for the OS (Windows or Linux) being used.

function compileBRILIA()
RootDir = findRoot;
if ispc
    OS = 'Win';
else
    OS = 'Linux';
end
TargetDir = fullfile(RootDir, 'Bin', OS);
[Success, Msg] = mkdir(TargetDir);
assert(Success > 0, '%s: Could not create directory "%s".\n  %s', mfilename, TargetDir, Msg)

compileMEX;

fprintf('%s: %s\n', mfilename, 'Compiling help texts.');
compileHelpText('Dir', RootDir, 'CheckSub', 'y', 'SaveTo', fullfile(RootDir, 'HelpText'), 'Overwrite', 'y');

fprintf('%s: %s\n', mfilename, 'Compiling BRILIA.m.');
mcc('-m', fullfile(RootDir, 'Src', 'BRILIA.m'), ...
    '-a', fullfile(RootDir, 'Src'), ...
    '-a', fullfile(RootDir, 'Tables'), ...
    '-a', fullfile(RootDir, 'HelpText'), ...
    '-a', fullfile(RootDir, 'Databases'));
if ~strcmp(TargetDir, cd)
    if ispc
        movefile('BRILIA.exe', TargetDir, 'f');
    elseif isunix
        movefile('BRILIA', TargetDir, 'f');
        movefile('run_BRILIA.sh', TargetDir, 'f');
    end
end

MatlabFiles = {'requiredMCRProducts.txt', 'mccExcludedFiles.log', 'readme.txt'};
for j = 1:length(MatlabFiles)
    if exist(MatlabFiles{j}, 'file')
        try
            delete(MatlabFiles{j});
        catch
        end
    end
end

zipBRILIA;