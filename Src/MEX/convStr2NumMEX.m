%convStr2NumMEX will convert a string to a number or a group of cell with 
%string to a group of cell with number/string
%
%  Num = convStr2NumMEX(Str)
%
%  INPUT
%    Str: a string of numbers separated by a non-decimal character
%         Note: Can accept scientific notations
% 
%  OUTPUT
%    Num: a scalar, vector, or []
%
%  EXAMPLE
%    Str = '[123 -123 1.23e2 -1.23e2 1.23e-2]';
%    Num = convStr2NumMEX(Str)
%    Num =
%      123.0000 -123.0000  123.0000 -123.0000    0.0123
%    
%    Num = convStr2NumMEX('34')
%    Num = 
%       34
%
%
%
