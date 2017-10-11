%zipBRILIA will zip up binary files in the BRILIA_ROOT_DIR/Bin/Win and
%BRILIA_ROOT_DIR/Bin/Linux folders.

function zipBRILIA
RootDir = findRoot;

LinuxDir = fullfile(RootDir, 'Bin', 'Linux');
try
    zip(fullfile(LinuxDir, 'BRILIA_Linux'), {'BRILIA', 'run_BRILIA.sh'}, LinuxDir);
catch
    warning('%s: Could not zip Linux files', mfilename);
end

WinDir = fullfile(RootDir, 'Bin', 'Win');
try
    zip(fullfile(WinDir, 'BRILIA_Win'), 'BRILIA.exe', WinDir);
catch
    warning('%s: Could not zip Windows files', mfilename);
end

