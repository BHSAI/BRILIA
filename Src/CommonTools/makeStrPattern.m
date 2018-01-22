%makeStrPattern is used to make a string like 'a,b,c,d,f' for use in repmat
%or any function required a string format.
%
%  Pattern = makeStrPattern(Str, Delimiter)  
%
%  INPUT
%    Str: cell array of char, or a char (each letter gets separated)
%    Delimiter: delimiter to use to separate the Str elements
%
%  OUTPUT
%    Pattern: a char with alternating Str elements and Delimiter
%
%  EXAMPLE
%    Str = {'a' 'b' 'c'};
%    Delimiter = ',';
%    Pattern = makeStrPattern(Str, Delimiter);
%    Pattern = makeStrPattern([Str{:}], Delimiter);
%    Pattern = 
%       'a,b,c'

function Pattern = makeStrPattern(Str, Delimiter)
if ~ischar(Delimiter)
    error('%s: Delimiter should be a char.')
end
if ischar(Str)
    Pattern = repmat(Delimiter, 1, length(Str)*2-1);
    Pattern(1:2:end) = Str;
elseif iscell(Str) && all(cellfun(@ischar, Str))
    Pattern = [sprintf(['%s' Delimiter], Str{1:end-1}) Str{end}];
else
    error('%s: Str should be a cell array of char or a char.')
end