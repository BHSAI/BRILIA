%chooseFromList will ask user to choose a value from a list N number of
%times. Returns the choice number and the choice. If it fails, will return
%an error.
%
%  [ChoiceNum, Choice] = chooseFromList(List)
%
%  [ChoiceNum, Choice] = chooseFromList(List, Param, Value, ...)hoice)
%
%  INPUT
%    List: Mx1 cell of string of a list of choices 
%
%     Param         Value (* Default)  Details
%     ------------  -----------------  --------------------------------
%     Attempt      * 1                 Max tries to select valid choice.
%     Default      * []                Return empty for failed choice.
%                  #                   Return this choice if user presses
%                                      Enter or there's no valid choice.
%     MultiSelect  * 'off'             Can only choose 1.
%                    'on'              Can choose multiple via '1:N' or
%                                        '1,2,3'
%  OUTPUT
%    ChoiceNum: integer value(s) for the List chosen
%    Choice: cell of values for the List chosen
%
%  EXAMPLE
%    List = {'What'; 'Is'; 'This'; 1322 ; 2332; 31423}
%    [ChoiceNum, Choice] = chooseFromList(List, 'Attempt', 5, 'Default', 1:3, 'Message', 'Choose from list: ')
%
function [ChoiceNum, Choice] = chooseFromList(varargin)
P = inputParser;
addRequired(P, 'List', @iscell)
addParameter(P, 'Attempt',     1,     @(x) isnumeric(x) && all(x >= 1) && all(mod(x, 1) == 0));
addParameter(P, 'Default',     [],    @(x) isempty(x) || (isnumeric(x) && all(x >= 1) && all(mod(x, 1) == 0)));
addParameter(P, 'MultiSelect', 'off', @(x) ischar(x) && ismember(lower(x), {'off', 'on'}));
addParameter(P, 'Message', 'Select #: ', @ischar);   
parse(P, varargin{:});
P = P.Results;

dispList(P.List);
ChoiceNum = -1;
j = 0;
while true
    UserChoice = input([P.Message ' '], 's');
    try
        UserNum = round(convStr2Num(UserChoice));
        if strcmpi(P.MultiSelect, 'off') && numel(UserNum) > 1
            fprintf('Cannot have multiple choices when MultiSelect is "off".\n');
        elseif isempty(UserNum) && ~isempty(P.Default)
            fprintf('Choosing default choice(s): %d\n', P.Default);
            ChoiceNum = P.Default;
            break
        elseif min(UserNum) >= 1 && max(UserNum) <= length(P.List)
            if strcmpi(P.MultiSelect, 'on')
                ChoiceNum = UserNum;
            else
                ChoiceNum = UserNum(1);
            end
            break
        end
    catch
    end
    fprintf('Invalid choice.\n')
    j = j + 1;
    if ChoiceNum(1) > 1; break; end
    if j >= P.Attempt
        fprintf('Failed to choose a valid choice in %d attempts.\n', j);
        if isempty(P.Default)
            ChoiceNum = [];
        else
            fprintf('Choosing default choice(s): %d\n', P.Default);
            ChoiceNum = P.Default;
        end
        break
    end
end
Choice = P.List(ChoiceNum);