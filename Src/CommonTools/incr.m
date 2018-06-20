%incr will increase a number by 1, but will set the max value to
%either a maximum number or a resetvalue.
%
%  Num = incr(Num)
%
%  Num = incr(Num, MaxNum)
%
%  Num = incr(Num, MaxNum, ResetValue)
%
%  INPUT
%    Num: scalar number
%    MaxNum: maximum number that it cannot exceed upon incrementing
%    ResetValue: the value to set Num to when the incrment is > MaxNum.
%      Default is the MaxNum.
%
%  OUTPUT
%    Num: Num+1 or 
%         MaxNum if ResetValue is empty or 
%         ResetValue if Num+1 > MaxNum


function Num = incr(Num, MaxNum, ResetValue)
Num = Num + 1;
if nargin >= 2 && Num > MaxNum
    if nargin == 3
        Num = ResetValue;
    else
        Num = MaxNum;
    end
end