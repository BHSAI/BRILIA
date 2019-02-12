%refactorVar will search through all *.m files in a directory and replace
%variable VAR1 to VAR2. This is not undoable, so please do save a backup
%before running this if needed. Before it refactors, it will show the user
%what variables will be refactored and this can be modified.
%
%  refactorVar(OldVar, NewVar)
%
%  refactorVar(OldVar, NewVar, Param, Value, ...)
%
%  INPUT
%    OldVar: old word to be replaced (case sensitive, whole word)
%    NewVar: new word that will replace OldVar
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
%  NOTE
%    The reason 'Exact' option is default to 'n' is because code variables
%    are recognized differently than word phrases. For instance, Exact =
%    'y' would be appropriate for replacing complex code phrase like
%    'Var(1:2)', but NOT 'en' - which would replace 'end' too. 
%
%  EXAMPLE
%    %To replace a word
%    refactorVar('this', 'that')
% 
%    %To replace a structure field
%    refactorVar('.fieldX', '.fieldY')
%
%    %To replace an exact phrase
%    refactorVar('Var(1:2)', 'Var(3:4)', 'Exact', 'y')
%
function refactorVar(OldVar, NewVar, varargin)
VarInfo = searchVar(OldVar, varargin{:});
if isempty(VarInfo)
    fprintf('No matches for %s\n', OldVar);
    return
end

%Display what lines are going to be replaced
while 1
    KeepLoc = ones(size(VarInfo, 1), 1, 'logical');
    DispInfo = cellfun(@(x, y, z) sprintf('%s | Line %d | %s', x, y, z), VarInfo(:, 2), VarInfo(:, 3), VarInfo(:, 4), 'un', 0);   
    dispList(DispInfo)
    DelInput = input('Select # to prevent making changes to these lines, or press ENTER to continue.\nWill ask for final confirmation. : ', 's');
    if isempty(DelInput); break; end
    KeepLoc(convStr2Num(DelInput)) = 0;
    VarInfo = VarInfo(KeepLoc, :);
    if isempty(VarInfo)
        fprintf('No replacement to be made for %s\n', OldVar);
        return
    end
end

%Final warning before refactoring
while 1
    ConfirmThis = input('WARNING: Refactoring variable is permanent. Continue? y or n: ', 's');
    if ~isempty(ConfirmThis)
        if strcmpi(ConfirmThis(1), 'y')
            break
        elseif strcmpi(ConfirmThis(1), 'n')
            fprintf('Aborted refactoring.\n');
            return
        end
    end
end

%Determine unique files to process
FullFileName = fullfile(VarInfo(:, 1), VarInfo(:, 2));
[~, ~, UnqNum] = unique(FullFileName);

%Proceed with refactoring
for f = 1:max(UnqNum)
    UnqIdx = find(UnqNum == f);
    FID = fopen(FullFileName{UnqIdx(1)}, 'r');
    if FID < 0
        warning('%s: Could not open "%s". Skipping search.', mfilename, FullFileName{UnqIdx(1)});
        continue
    end
    TXT = textscan(FID, '%s', 'Delimiter', '\n', 'WhiteSpace', '');
    fclose(FID);

    LineIdx = cell2mat(VarInfo(UnqIdx, 3));
    TXT{1}(LineIdx) = replaceCodeVar(TXT{1}(LineIdx), OldVar, NewVar, varargin{:});
    
    FID = fopen(FullFileName{UnqIdx(1)}, 'w');
    fprintf(FID,'%s\n', TXT{1}{:});
    fclose(FID);
end