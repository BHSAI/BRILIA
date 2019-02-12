%replaceCodeVar will search for an exact word and case-senstive match to a
%code variable. This differs from strfind and regexp in that it yields the
%start index of a code variable while accounting for the structure "."
%notations to enable the replacement of structure names as well.
%
%  Idx = replaceCodeVar(Txt, OldVar, NewVar)
%
%  INPUT
%    Txt: 1xN char line of a code, or a Mx1 cell of 1xN char line of codes
%    OldVar: old code variable.
%    NewVar: new code variable to replace OldVar.
%
%  OUTPUT
%    IdxMat: Qx1 index of 1st letter of Var in Q occurrences of a line of
%      code, or a Mx1 cell of Qx1 indices
%
%  NOTE
%    A '.' is needed before a structure field name to replace that field.  
%
%  EXAMPLE
%    Txt{1, 1} = '   A.var1 = 2';
%    Txt{2, 1} = '   Test.A = 3';
%    Txt{3, 1} = '   A.var1 = (A.var1+3)+A.var1/A.var1';
%
%    %To replace a normal variable or code word
%    NewTxt = replaceCodeVar(Txt, 'A', 'B');
%    NexTxt = 
%      3x1 cell
%        '   B.var1 = 2'
%        '   Test.A = 3'
%        '   B.var1 = (B.var1+3)+B.var1/B.var1'
%
%    %To replace a structure field name with another field name
%    NewTxt = replaceCodeVar(NewTxt, '.A', '.C'); 
%      3x1 cell
%        '   B.var1 = 2'
%        '   Test.C = 3'
%        '   B.var1 = (B.var1+3)+B.var1/B.var1'
%
function Txt = replaceCodeVar(Txt, OldVar, NewVar, varargin)
if iscell(Txt)
    Txt = cellfun(@(x) replaceCodeVarPerLine(x, OldVar, NewVar, varargin{:}), Txt, 'un', 0);
else
    Txt = replaceCodeVarPerLine(Txt, OldVar, NewVar, varargin{:});
end

function NewTxt = replaceCodeVarPerLine(Txt, OldVar, NewVar, varargin)
VarIdx = findCodeVar(Txt, OldVar, varargin{:});
if isempty(VarIdx) 
    NewTxt = Txt;
    return
end

if VarIdx(1) > 1
    InitTxt = Txt(1:VarIdx(1)-1);
else
    InitTxt = '';
end
NewTxt = cell(1, numel(VarIdx));
for j = 1:length(VarIdx)-1
    NewTxt{j} = [NewVar Txt(VarIdx(j)+length(OldVar):VarIdx(j+1)-1)];
end
NewTxt{end} = [NewVar Txt(VarIdx(end)+length(OldVar):end)];
NewTxt = horzcat(InitTxt, NewTxt{:});