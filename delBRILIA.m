%delAllPaths will delete paths assocatiated with BRILIA.
function delBRILIA
if ~isempty(which('findRoot'))
    RootDir = findRoot;
    warning('off','all')
    rmpath(genpath(RootDir));
    warning('on','all')
end