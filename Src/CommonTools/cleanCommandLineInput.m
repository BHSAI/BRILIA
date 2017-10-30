%cleanCommandLineInput will remove the '-' in in a cell array of string inputs
%such as varargin, and convert numerical strings to numbers. This is used
%mainly for command line usage of a function where a '-' is used in front
%of the parameter name.
%
%  INPUT
%    CellStr: a cell array of param/value pairs like from varargin
%
%  OUTPUT
%    CellStr: CellStr with the leading '-' removed from characters that are
%      not negative numbers
%
%  EXAMPLE
%    varargin = {'-in', 'in.txt', '--out', 'out.txt', '-count', '-1'};
%    varargin = cleanCommandLineInput(varargin{:})
%    varargin = 
%       'in'    'in.txt'    'out'    'out.txt'    'count'    [-1]
function CellStr = cleanCommandLineInput(varargin)
if length(varargin) == 1 && iscell(varargin{1})
    CellStr = varargin{1};
else
    CellStr = varargin;
end

CharLoc = find(cellfun(@ischar, CellStr));
for k = 1:length(CharLoc)
    %Check if this is a word with leading dashes
    InputStr = CellStr{CharLoc(k)};
    NonDashLoc = regexp(InputStr, '[^-]', 'once');
    if isempty(sscanf(InputStr(NonDashLoc), '%f'))
        CellStr{CharLoc(k)} = InputStr(NonDashLoc:end);
        continue
    end
    
    %Check if this a string that can be a number
    if isempty(regexpi(InputStr, '[^0-9\]\[\-\.\:\,]'))
        try
            NumVal = convStr2Num(InputStr);
            if ~isempty(NumVal)
                CellStr{CharLoc(k)} = NumVal;
                continue;
            end
        catch
        end
    end
end