%compileMex gets rid of the nuances of trying to add on the include files
%based on the custom header files included in the main Mex function . 
%
%  compileMex(FileName)
%
%  compileMex(FileName, varargin)
%
%  compileMex(..., 'genmfile')
%
%  INPUT
%     FileName: FileName you want to compile
%     varargin: 1xM cell of any inputs wanted for the mex compiler. Note
%     that to specify a include folder, use the "-I<Path>" notation of the
%     mex compiler.
%     '-gennmfile': will generate .m files based on the comments ABOVE the
%        first #include statement. The .m files will be at the same place
%        as the .mex files.
%
%  OUTPUT
%    Compiled mex code. Same as doing:
%      >> mex FileName -I<srcfile.cpp> -I<srcfile.cpp> 
%
%  NOTE
%    This does NOT look for dependencies of the source file header files.
%    It only does it for the main MEX function. The source file should only
%    have #include "mex.h" and/or <*.h> files.
%
%    For .hpp files, will look for .cpp files.
%    For .h files, will look for .c files.
%    
%   compileMex('convAA2PropMEX.cpp')
function compileMex(FileName, varargin)
GenMFileLoc = contains(varargin, 'genmfile', 'ignorecase', true); 
DoGenMFile = any(GenMFileLoc);
if DoGenMFile
    varargin(GenMFileLoc) = [];
end

OutIdx = find(contains(varargin, 'outdir', 'ignorecase', true), 1);
if ~isempty(OutIdx)
    OutDir = varargin{OutIdx + 1};
    varargin(OutIdx:OutIdx+1) = [];
else
    OutDir = [];
end    

IncLoc = startsWith(varargin, '-I');
SrcPath = cellfun(@(x) x(3:end), strtrim(varargin(IncLoc)), 'un', false);

FileName = which(FileName);
[FP, FN] = parseFileName(FileName);
assert(~isempty(FP), '%s: Could not find the file "%s".', mfilename, FileName);
FullName = fullfile(FP, FN);

if isempty(OutDir); OutDir = fileparts(FullName); end
InclFiles = findIncludeFiles(FullName, SrcPath);
UnqSrcPath = unique(cellfun(@(x) ['-I' fileparts(x)], InclFiles, 'un', false)');
mex(FullName, InclFiles{:}, UnqSrcPath{:}, '-outdir', OutDir, varargin{:});

if DoGenMFile
    genMexMFile(FullName);
end

%findIncludeFiles will search for all #include source files that is
%associated with the input file given as input. It will look either in a
%specified SrcPath or any subdir where the input file is. 
function [IncFiles, SrcPath] = findIncludeFiles(FullName, SrcPath, varargin)
%Find a SrcPath to look for source codes
FilePath = parseFileName(FullName);
if nargin == 1 || isempty(SrcPath)
    SrcPath = fullfile(FilePath, '**', filesep); %default to any sub dir
end
if ischar(SrcPath)
    SrcPath = strsplit(SrcPath, ';');
end
SrcPath = unique(strtrim(SrcPath));
AddFiles = extractIncludeFiles(FullName, SrcPath);
IncFiles = AddFiles;
while ~isempty(AddFiles)
    TmpFiles = cell(1, length(AddFiles));
    for j = 1:length(TmpFiles)
        TmpFiles{j} = extractIncludeFiles(AddFiles{j}, SrcPath);
    end
    TmpFiles = unique([TmpFiles{:}]);
    IncFiles = unique([IncFiles TmpFiles]);
    AddFiles = setdiff(TmpFiles, IncFiles);
end

%extractIncludeFiles will extract the C/H source codes that are referred to
%by the #include statements in the first C source code. This will NOT
%recursively search for all source codes.
function SrcFiles = extractIncludeFiles(FullName, SrcPath)
%Extract all the #include statements
[FID,  MSG] = fopen(FullName, 'r');
assert(FID > 0, '%s: Could not open file "%s". \n  MSG: %s.', mfilename, MSG);
TXT = textscan(FID, '%s', 'delimiter', '\n');
fclose(FID);
IncludeLine = TXT{1}(cellfun(@(x) ~contains(strtrim(x), {'<', 'mex.h', 'cytpe.h', '//#include'}) && contains(strtrim(x), '#include'), TXT{1}));

%Get just the include source code files
SrcFiles = cell(1, length(IncludeLine));
for j = 1:length(SrcFiles)
    QuoteIdx = find(IncludeLine{j} == '"');
    assert(numel(QuoteIdx) == 2, '%s: Invalid syntax for the #include "name.hpp" statement. Recheck "%s".', mfilename, IncludeLine{j});
    IncludeName = IncludeLine{j}(QuoteIdx(1)+1:QuoteIdx(2)-1);
    DotLoc = find(IncludeName == '.', 1, 'last');
    IncludeName(DotLoc+1) = 'c';
    SrcFiles{j} = IncludeName;
end

%Get the full path to the source code files
for j = 1:length(SrcFiles)
    for k = 1:length(SrcPath)
        Src = dir(fullfile(SrcPath{k}, SrcFiles{j}));
        if ~isempty(Src); break; end
    end
    if isempty(Src)
        warning('%s: Could not find the file "%s".', fullfile(SrcFiles{j}));
        SrcFiles{j} = [];
    else
        SrcFiles{j} = fullfile(Src.folder, Src.name);
    end
end
SrcFiles = unique(SrcFiles(~cellfun('isempty', SrcFiles)));
