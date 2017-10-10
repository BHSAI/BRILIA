%DUPLICATE of convStr2Num to remove dependencies.
%
%convStr2Num will convert a string to a number, which can handle integer
%range notations (ex: '1:4,-3,20'). 
%
%  Num = convStr2Num(Str)
%
%  Num = convStr2Num(Str, 'unique')
%
%  INPUT
%    Str: a numeric or range string
%    'unique': unique and sorts the nubmer in ascending order
%  
%  OUTPUT
%    Num: a scalar or vector of numbers specfieid by the string
%         OR, [] if it fails at conversion.
%
%  EXAMPLE
%    Str = '1:4,-3,3:5';
%    Num = convStr2Num(Str)
%    Num =
%         1   2   3   4  -3   3   4   5
%    
%    Num = convStr2Num(Str, 'unique')
%    Num =
%        -3   1   2   3   4   5

function Num = convStr2NumT(Str, Option)
GetUnq = false;
if nargin == 2 && strcmpi(Option, 'unique')
    GetUnq = true;
end

Vec = sscanf([Str ','],'%f%c');
if isempty(Vec)
    Num = [];
    return
elseif numel(Vec) == 1
    Num = Vec;
    return
end
Num = [];
s1 = 1;
while s1 <= numel(Vec)
    s2 = s1+2*(Vec(s1+1)==58); % char(58) is ':'
    Num = cat(2, Num, Vec(s1):Vec(s2));
    s1 = s2+2;
end

if GetUnq
    Num = unique(Num);
end