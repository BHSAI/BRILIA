%addBRILIA will add all folders and subfolders associated with BRILIA.
function addBRILIA
addpath(genpath(fileparts(mfilename('fullpath'))));
rmpath(genpath(fullfile(findRoot, 'Bin')));
