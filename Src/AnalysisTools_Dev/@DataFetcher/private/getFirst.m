%getFirst will get the first name of a "|"-separated list of char items.
%
%  FirstOnly = getFirst(Name)
%
%  INPUT
%    Name: char or cell array of chars
%
%  OUTPUT
%    FirstOnly: char or cell array of chars of the first item of a list
%
%  EXAMPLE
%    Name = 'A|B|C|D';
%    FirstOnly = getFirst(Name);
%    FirstOnly =
%      'A'
%
function FirstOnly = getFirst(Name)
if iscell(Name)
    FirstOnly = cell(size(Name));
    parfor j = 1:length(Name)
        FirstOnly{j} = getFirst_Calc(Name{j});
    end
elseif ischar(Name)
    FirstOnly = getFirst_Calc(Name);
else
    error('%s: Name must be a char or a cell array of char.', mfilename);
end

function FirstOnly = getFirst_Calc(Name)
DivLoc = find(Name == '|', 1);
if isempty(DivLoc)
    FirstOnly = Name;
else
    FirstOnly = Name(1:DivLoc-1);
end

