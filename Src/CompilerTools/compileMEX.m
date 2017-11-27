%compileMEX will compile MEX functions of BRILIA, and store them in the
%correct folders.
function compileMEX
RootDir = findRoot;
if ispc
    OS = 'Win';
else
    OS = 'Linux';
end

MexDir = fullfile(RootDir, 'Src', 'MEX');
MexTargetDir = fullfile(RootDir, 'Src', 'MEX', OS);
[Success, Msg] = mkdir(MexTargetDir);
assert(Success > 0, '%s: Could not create directory "%s".\n  %s', mfilename, MexTargetDir, Msg)
MexFiles = dir(fullfile(MexDir, '*.cpp'));
for f = 1:length(MexFiles)
    mex(fullfile(MexDir, MexFiles(f).name), '-outdir', MexTargetDir);
end
