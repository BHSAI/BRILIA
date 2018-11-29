%zipBRILIA will zip up binary files in the BRILIA_ROOT_DIR/Bin/Win and
%BRILIA_ROOT_DIR/Bin/Linux folders.

function zipBRILIA(varargin)
DoAll = false;
if nargin >= 1 && islogical(varargin{1})
    DoAll = varargin{1};
end

RootDir = findRoot;
if ispc || DoAll
    WinDir = fullfile(RootDir, 'Bin', 'Win');
    try
        zip(fullfile(WinDir, sprintf('BRILIAv%s_Win', BRILIA('version'))), 'BRILIA.exe', WinDir);
    catch
        warning('%s: Could not zip Windows files', mfilename);
    end

elseif isunix || DoAll
    LinuxDir = fullfile(RootDir, 'Bin', 'Linux');
    try
        zip(fullfile(LinuxDir, sprintf('BRILIAv%s_Linux', BRILIA('version'))), {'BRILIA', 'run_BRILIA.sh'}, LinuxDir);
    catch
        warning('%s: Could not zip Linux files', mfilename);
    end
end

