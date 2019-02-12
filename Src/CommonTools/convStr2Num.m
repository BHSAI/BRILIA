%convStr2Num will convert a string to a number, which can handle integer
%range notations (ex: '1:4,-3,20'). 
%
%  Num = convStr2Num(Str)
%
%  INPUT
%    Str: a numeric or range string
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
%  NOTE
%    This does not perform different increments other than 1.
%    Str = '1:-1:4' % will NOT work!

function Num = convStr2Num(Str)
if isempty(Str)
    Num = [];
    return
end

Str = strrep(Str, ',', ' '); %better to shift , to space as sscanf will skip spaces
t = (Str(1) == '[');
Vec = sscanf([Str(1+t:end), ','],'%f%c');
if numel(Vec) <= 1
    Num = Vec;
    return
end

%CODING_NOTE: Growing matrices is still faster most of the time since it's small
Num = [];
a = 1;
while a <= numel(Vec)
    t = (Vec(a+1) == 58);
    b = a + 2*t; % char(58) is ':'
    Num = cat(2, Num, Vec(a):Vec(b));
    a = b + 2;
end