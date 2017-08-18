%regexp2 extends regexp basic capabilities of regexp, but will return
%overlapping matches, which is shown example below:
%
%  Str = 'AAAAA'
%  Expr = 'AA'
%  StartLoc = regexp(Str, Expr) %Does NOT find overlapping AA patterns
%  StartLoc = 
%        1  3
%
%  StartLoc = regexp2(Str, Expr) 
%  StartLoc =
%        1  2  3  4

function varargout = regexp2(Str, Expr)
StrIdx = regexp(Str, Expr, 'once');
MatchIdx = zeros(1, length(Str), 'logical');
t = 0;
while ~isempty(StrIdx)
    t = t + StrIdx;
    MatchIdx(t) = 1;
    StrIdx = regexp(Str(t+1:end), Expr, 'once');
end
varargout{1} = find(MatchIdx);
