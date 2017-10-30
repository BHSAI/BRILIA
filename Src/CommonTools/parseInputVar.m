%parseInputVar extracts a parameter of a param-value pair within varargin.
%This avoids having to use inputParser (P), and extracting the value in a
%lengthy structure.
%
%  ParamVal = parseInputVar(ParamName, Default, Validator, varargin{:});
%
%  INPUT
%    ParamName: Name of parameter
%    Default: Default parameter value
%    Validator: 1 or a function handle (ex: @(x)islogical(x) ) for ensuring
%      the value is valid
%    varargin: the varargin input cell of another function
%
%  OUTPUT
%    ParamVal: Value of the parameter
%
%  NOTE
%    This works only for varargin with param-value pairs. If there are
%    optional or special inputs, those must be parsed separately and
%    removed prior to passing into this function.
%
%  EXAMPLE (parseInputVar - this code)
%    varargin = {'a', 1, 'b', 2};
%    P = struct;
%    [a, P] = parseInputVar('a', 0, @(x)isnumeric(x), varargin{:}, P);
%    [c, P] = parseInputVar('c', 3, @(x)isnumeric(x), varargin{:}, P);
%
%  EXAMPLE (inputParser - messy to deal with)
%    varargin = {'a', 1, 'b', 2};
%    P = inputParser;
%    P.KeepUnmatched = 1;
%    P.addParameter('a', 0, @(x) isnumeric(x));
%    P.addParameter('c', 3, @(x) isnumeric(x));
%    parse(P, varargin{:});
%    a = P.Results.a;
%    c = P.Results.c;

function [ParamVal, P] = parseInputVar(ParamName, Default, Validator, varargin)
%Determine varargin before any input structure is given at the end
if ~isempty(varargin) && isstruct(varargin{end})
    EndLoc = length(varargin) - 1;
else
    EndLoc = length(varargin);
end

%Extract the value for this paramater
VarLoc = 0;
for j = 1:2:EndLoc
    if strcmpi(varargin{j}, ParamName)
        if j+1 > EndLoc
            error('%s: Could not find value for (%s) - exceeded length(varargin).', mfilename, ParamName);
        else
            VarLoc = j+1;
            ParamVal = varargin{VarLoc};
        end
        break;
    end
end
if VarLoc == 0
    ParamVal = Default;
end

%Validate input type
if isa(Validator, 'function_handle')
    BadFuncStr = {'delete', 'copyfile', 'movefile', 'mkdir', 'rmdir'};
    ValidatorStr = func2str(Validator);
    ValidatorStrOnly = regexpi(ValidatorStr, '[a-z0-9]+', 'match');
    if any(ismember(BadFuncStr, ValidatorStrOnly))
        error('%s: Validator "%s" cannot be a system function handles.', mfilename, ValidatorStrOnly.name); 
    end
    try
        ValidatorResult = Validator(ParamVal);
    catch
        error('%s: Validator "%s" could not be evaluated.', mfilename, ValidatorStr);
    end
    if numel(ValidatorResult) ~= 1 || ~islogical(ValidatorResult)
        error('%s: Validator "%s" must return one true/false result.', mfilename, ValidatorStr);
    end
    if ~ValidatorResult
        error('%s: Input parameter "%s" failed validation by "%s".', mfilename, ParamName, func2str(Validator));
    end
else
    error('%s: Validator must be a function handle.', mfilename);
end

%Modify input structure if given and if it's an output
if nargout > 1
    if isstruct(varargin{end})
        P = varargin{end};
        P.(ParamName) = ParamVal;
    else
        P = struct('ParamName', ParamVal);
    end
end