%compileMEX will compile MEX functions of BRILIA, and store them in the
%correct folders.
function compileMEX
RootDir = findRoot;
if ispc
    OS = 'Win';
else
    OS = 'Linux';
end

%Create the target location
MexTargetDir = fullfile(RootDir, 'Src', 'MEX', OS);
[Success, Msg] = mkdir(MexTargetDir);

%Find all mex cpp files under the root dir
MexFiles = dir(fullfile(RootDir, '**', '*.cpp'));
MexFileNames = arrayfun(@(x) fullfile(x.folder, x.name), MexFiles, 'unif', false);
assert(Success > 0, '%s: Could not create directory "%s".\n  %s', mfilename, MexTargetDir, Msg)
for f = 1:length(MexFileNames)
    mex(MexFileNames{f}, '-outdir', MexTargetDir);
end
