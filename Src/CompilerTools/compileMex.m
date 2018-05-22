%compileMex will compile MEX functions of BRILIA, and store them in the
%correct folders. This will also generate the associated .m files used for
%Matlab help functions.
%
%  compileMex
%
%  compileMex(FileNames)
%
%  INPUT
%    FileNames: file name string or cell array of strings. If empty, will
%      do all mex files within BRILIA
function compileMex(FileNames)
if nargin == 0
    RootDir = findRoot;
    Files = dir(fullfile(RootDir, '**', '*.cpp'));
    FileNames = arrayfun(@(x) fullfile(x.folder, x.name), Files, 'unif', false);
else
    if ischar(FileNames)
        FileNames = {FileNames};
    end
    for f = 1:length(FileNames)
        [FilePath, FileName] = parseFileName(FileNames{f});
        FileNames{f} = fullfile(FilePath, FileName);
    end
end

for f = 1:length(FileNames)
    OutDir = fileparts(FileNames{f});
    mex(FileNames{f}, '-outdir', OutDir);
    compileMexMFiles(FileNames{f})
end
