%decr will decrease a number by 1, but will set the min value to
%either a minimum number or a resetvalue.
%
%  Num = decr(Num)
%
%  Num = decr(Num, MinNum)
%
%  Num = decr(Num, MinNum, ResetValue)
%
%  INPUT
%    Num: scalar number
%    MinNum: minimum number that it cannot exceed upon decrementing
%    ResetValue: the value to set Num to when the decrement is < MinNum.
%      Default is the MinNum.
%
%  OUTPUT
%    Num: Num-1 or 
%         MinNum if ResetValue is empty or 
%         ResetValue if Num-1 < MinNum

function Num = decr(Num, MinNum, ResetValue)
Num = Num - 1;
if nargin >= 2 && Num < MinNum
    if nargin == 3
        Num = ResetValue;
    else
        Num = MinNum;
    end
end