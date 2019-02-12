%searchVar will search through all *.m files in a directory for a variable
%and save information that helps locate the occurrence of this variable.
%
%  Result = searchVar(Var)
%
%  Result = searchVar(Var, Param, Value, ...)
%
%  INPUT
%    Var: string to search within the m files. Only full word matches are
%      allowed, and it is case-senstive.
%
%    Param           Value          Description
%    ---------       ------------   ---------------------------
%    'Dir'           [path]         Folder name to check for m files. 
%                                     If empty, will use current directory.
%    'CheckSub'      'y' 'n'        Yes/No to checking subdirectories too
%                                     Default is 'y'
%    'Exact'         'n' 'y'        No/Yes to exact matching of phrases.
%                                     Default is 'n' since code variable is
%                                     not the same as an exact phrase. See
%                                     Note. 
%
%  OUTPUT
%    Result: a Mx4 cell matrix containing the following information:
%      Col1 - m file path
%      Col2 - m file name
%      Col3 - line number
%      Col4 - code text
%      
%  NOTE
%    The reason 'Exact' option is default to 'n' is because code variables
%    are recognized differently than word phrases. For instance, Exact =
%    'y' would be appropriate for replacing complex code phrase like
%    'Var(1:2)', but NOT 'en' - which would replace 'end' too. 
%
%  EXAMPLE
%    Result = searchVar('makeString', 'Dir', cd, 'CheckSub', 'n')

function Result = searchVar(Var, varargin)
P = inputParser;
P.CaseSensitive = 0;
P.PartialMatching = 1;
addParameter(P, 'Dir', pwd, @ischar);
addParameter(P, 'CheckSubDir', 'y', @(x) ischar(x) && startsWith(x, {'y', 'n'}, 'ignorecase', true));
addParameter(P, 'Exact', 'n', @(x) ischar(x) && startsWith(x, {'y', 'n'}, 'ignorecase', true));
parse(P, varargin{:});
Dir = P.Results.Dir;
CheckSubDir = P.Results.CheckSubDir;
Exact = P.Results.Exact;

%Determine the subfolders too
if startsWith(CheckSubDir, 'y', 'ignorecase', true)
    MDir = dir(fullfile(Dir, '**', '*.m'));
else
    MDir = dir(fullfile(Dir, '*.m'));
end

Result = cell(length(MDir), 1);
for f = 1:length(MDir)
    FID = fopen(fullfile(MDir(f).folder, MDir(f).name), 'r');
    if FID < 0
        warning('%s: Could not open "%s". Skipping search.', mfilename, MDir(f).name);
        continue
    end
    TXT = textscan(FID, '%s', 'Delimiter', '\n', 'WhiteSpace', '');
    fclose(FID);
    LineIdx = find(~cellfun(@isempty, findCodeVar(TXT{1}, Var, Exact)));
    Result{f} = cell(length(LineIdx), 4);
    Result{f}(:, 1) = {MDir(f).folder};
    Result{f}(:, 2) = {MDir(f).name};
    Result{f}(:, 3) = num2cell(LineIdx);
    Result{f}(:, 4) = TXT{1}(LineIdx);
end
Result = vertcat(Result{:});