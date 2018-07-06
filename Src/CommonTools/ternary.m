%ternary perform ternary operation for cases where one does if this is
%true, return value 1, otherwise return value 2. Ternary operators are
%often used in other lanugages like c++ and java, such as
%value = (k < 0) ? 1 : 2;
%
%  Value = ternary(Equality, TrueVal, FalseVal);
%
%  Value = ternary(FuncHandle, TrueVal, FalseVal, varargin);
%
%  INPUT
%    Equality: logical 1 or 0 from a equality test. 
%    FuncHandle: function handle to compute the equality. Use the 4+ inputs
%       to supply the variable to the function handle.
%    varargin: inputs to the FuncHandle equality.
%    
%  OUTPUT
%    Value: either TrueVal or FalseVal
%
%  EXAMPLE
%    x = 2;
%    Value = ternary(x < 2, 1, 3)
%    Value =
%             3
%
%    Value = ternary(@(x) isempty(x) || x < 2, 1, 3, [])
%    Value =
%             1
%  
%
function Value = ternary(Equality, Value, FalseVal, varargin)
if isa(Equality, 'function_handle')
    Test = Equality(varargin{:});
elseif islogical(Equality)
    Test = Equality;
else
    error('%s: Input must be a logical value or function handle that returns a logical value.', mfilename);
end
if ~Test
    Value = FalseVal;
end