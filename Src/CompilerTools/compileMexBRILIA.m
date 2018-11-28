%compileMexBRILIA will compile MEX functions of BRILIA, and store them in
%the correct folders. This will also generate the associated .m files used
%for Matlab help functions.
%
%  compileMexBRILIA
%
function compileMexBRILIA
RootDir = findRoot;
MexFiles = dir(fullfile(RootDir, '**', '*MEX.cpp'));
MexFiles = fullfile({MexFiles.folder}, {MexFiles.name});
clear mex;  %#ok<CLMEX>  %Just include as otherwise could lead to write permission on compile
for f = 1:length(MexFiles)
    compileMex(MexFiles{f}, 'makeMFile', '-outdir', fileparts(MexFiles{f}));
end