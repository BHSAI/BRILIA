%cleanCommandLineInput will clean the variable inputs given to a function
%from a matlab or command line environment. From the command line, inputs
%are considered char until translated to a number. This code does these:
%1) checks for equal number of quotes ". Useful for filenames with spaces.
%2) removes any '-' that precedes a parmaeter option (eg '-species'), but
%   not a '-' preceding a negative number (eg '-1.23')
%3) splits the inputs into individual input cells
%4) converts "" to ''
%5) converts strings to numbers 
%
%  INPUT
%    CellStr: a cell array of param/value pairs like from varargin
%
%  OUTPUT
%    CellStr: cleaned and translated CellStr ready for use in matlab
%
%  EXAMPLE
%    varargin = {'-in' '"my fld\in.txt"' '--out' 'out.txt' '-ct' '-1:3'}
%    varargin = cleanCommandLineInput(varargin{:})
%    varargin =
%        'in'  '"my fld\in.txt"'  'out'  'out.txt'  'ct'  [1×5 double]
%
%    Input = '-in  "my folder\in.txt"  --out  out.txt  -ct  -1:3'
%    varargin = cleanCommandLineInput(Input)
%    varargin =
%        'in'  '"my fld\in.txt"'  'out'  'out.txt'  'ct'  [1×5 double]
%
function CellStr = cleanCommandLineInput(varargin)
if length(varargin) == 1 && ischar(varargin{1}) %input was a unparsed string
    StrInput = varargin{1};
    QuoteIdx = find(StrInput == '"');
    if mod(QuoteIdx, 2) > 0
        CellStr = {};
        warning('%s: Unequal number of quotes.', mfilename);
        return
    end
    for q = 1:2:length(QuoteIdx)
        StrInput(QuoteIdx(q):QuoteIdx(q+1)) = strrep(StrInput(QuoteIdx(q):QuoteIdx(q+1)), ' ', '^'); %Workaround, convert space to ^ within quotes to prevent split next
    end
    CellStr = strrep(strrep(regexp(StrInput, '\s+', 'split'), '^', ' '), '"', ''); %split by space, restore ^ to space, remove quotes
else
    CellStr = varargin;
    CharLoc = cellfun('isclass', CellStr, 'char');
    CellStr(CharLoc) = strrep(CellStr(CharLoc), '"', '');
end

%Determine if this is a matlab-style input, like 'function_name(a, b, c)'
if numel(CellStr) == 1 && ischar(CellStr{1}) && contains(CellStr{1}, {'(', ')'})
    StrOnly = regexprep(CellStr{1}, {' ', ')', ''''}, ''); %Remove spaces
    CellStr = strsplit(StrOnly, {'(', ')', ','});
end

CharIdx = find(cellfun('isclass', CellStr, 'char'));
for k = 1:length(CharIdx)
    Str = CellStr{CharIdx(k)};
    NonDashIdx = regexp(Str, '[^-]', 'once');
    if isempty(sscanf(Str(NonDashIdx), '%f'))
        CellStr{CharIdx(k)} = Str(NonDashIdx:end);
        continue
    end
    
    %Check if this a string that can be a number
    if isempty(regexpi(Str, '[^0-9\]\[\-\.\:\,]'))
        try
            NumVal = convStr2Num(Str); %Do not use the MEX version, as it will not convert 1:3 into 1,2,3.
            if ~isempty(NumVal)
                CellStr{CharIdx(k)} = NumVal;
                continue
            end
        catch
        end
    end
end