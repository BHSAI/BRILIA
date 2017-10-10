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
%    a = parseInputVar('a', 0, @(x)isnumeric(x), varargin{:});
%    c = parseInputVar('c', 3, @(x)isnumeric(x), varargin{:});
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

function ParamVal = parseInputVar(ParamName, Default, Validator, varargin)
VarLoc = strcmpi(varargin, ParamName);
if any(VarLoc)
    VarIdx = find(VarLoc) + 1;
    if VarIdx(1) > length(varargin)
        error('%s: Could not find value for (%s) - exceeded length(varargin).', mfilename, ParamName);
    end
    ParamVal = varargin{VarIdx(1)};
else
    ParamVal = Default;
end

if isa(Validator, 'function_handle')
    BadFuncStr = 'delete|copyfile|mkdir|rmdir';
    if ~isempty(regexpi(func2str(Validator), BadFuncStr))
        error('%s: Cannot accept function handles with these words: %s.', mfilename, BadFuncStr); 
    end
    ValidatorResult = Validator(ParamVal);
elseif Validator == 1
    ValidatorResult = 1;
else
    error('%s: Validator is invalid', mfilename);
end

if ~islogical(ValidatorResult)
    error('%s: Validator must return a single logical value.', mfilename);
elseif ~min(ValidatorResult)
    error('%s: Input parameter (%s) failed validation: %s.', mfilename, ParamName, func2str(Validator));
end