%findCodeVar will search for an exact word and case-senstive match to a
%code variable. This differs from strfind and regexp in that it yields the
%start and end index of a code variable, and accounts for the structure "."
%notations to enable the replacement of structure names as well.
%
%  Idx = findCodeVar(Txt, Var, Exact)
%
%  INPUT
%    Txt: 1xN char line of a code, or a Mx1 cell of 1xN char line of codes
%    Var: code variable to search for
%    Exact ['y' 'n']: yes or no to exact phrase matching
%
%  OUTPUT
%    IdxMat: Qx1 index of 1st letter of Var in Q occurrences of a line of
%      code, or a Mx1 cell of Qx1 indices
%
%  EXAMPLE
%    Txt{1, 1} = '   A.var1 = 2';
%    Txt{2, 1} = '   Test.A = 3';
%    Txt{3, 1} = '   A.var1 = (A.var1+3)+A.var1/A.var1';
%
%    %To find a normal variable or code word
%    VarIdx = findCodeVar(Txt, 'A');
%    VarIdx = 
%      3x1 cell array
%        [4]
%        []
%       [4 14 24 31]
%
%    %To find a structure field
%    VarIdx = findCodeVar(Txt, '.A');
%    VarIdx = 
%      []
%      [8]
%      []
%
function VarIdx = findCodeVar(Txt, Var, Exact)
if nargin < 3
    Exact = 'n'; 
end
if iscell(Txt)
    VarIdx = cellfun(@(x) findCodeVarPerLine(x, Var, Exact), Txt, 'un', 0);
else
    VarIdx = findCodeVarPerLine(Txt, Var, Exact);
end

function VarIdx = findCodeVarPerLine(TxtLine, Var, Exact)
VarIdx = strfind(TxtLine, Var);
if strcmpi(Exact(1), 'y'); return; end
if Var(1) == '.' %Fixing a structure field name
    LeftExclPat = '\W'; %Left of '.' cannot be anything other than a-zA-Z0-9
else
    LeftExclPat = '[^\s\W]|^\.'; %Left of Var cannot be anything other than a space or non-word symbol and '.'
end
for j = 1:length(VarIdx)
    %Look for invalid left of code var
    if VarIdx(j) > 1
        if ~isempty(regexp(TxtLine(VarIdx(j)-1), LeftExclPat, 'once'))
            VarIdx(j) = 0;
        end
    end
    
    %Look for invalid right of code var. "Var." is valid for struct var
    if VarIdx(j)+length(Var) <= length(TxtLine)
        if ~isempty(regexp(TxtLine(VarIdx(j)+length(Var)), '[^\.\s\W]', 'once'))
            VarIdx(j) = 0;
        end
    end
end
VarIdx = VarIdx(VarIdx > 0);