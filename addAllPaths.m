%addAllPaths is a simple 1-line code for adding all of folders and
%subfolders to the matlab path. 
%
%  NOTE
%    This code will be removed in future release. It's kept now since
%    it  helps with debugging and installing BRILIA paths quickly.
function addAllPaths
if ~isdeployed
    addpath(genpath(cd));
end

