%cmprFunctions will compare two functions given the same inputs to ensure
%they return the same outputs.
%
%  Summary = cmprFunctions(Func1, Func2, Nout, Iter, Inputs...)
%
%  INPUT
%    Func1: 1st function name or function handle
%    Func2: 2nd function name of function handle
%    Nout:  number of output arguments
%    Iter:  number of itereations to test
%    Inputs...: any comma-separated list of inputs for Func1 and Func2 that
%      will get passed on as Func1(varargin{:})
%
%  OUTPUT
%    Summary: structure storing the results of the function comparisons
%      .Time1: time to run Func1 Iter number of times
%      .Time2: time to run Func2 Iter number of times
%      .AssertTest: 1xNout logical array. 1 = match output, 0 = mismatch
%
%  EXAMPLE
%    Summary = cmprFunctions('convStr2NumMEX', 'convStr2Num', 1, 1000, '[10, 20, 234.4, -323.123]')

function Summary = cmprFunctions(Func1, Func2, Nout, Iter, varargin)

if ischar(Func1)
    Func1 = str2func(Func1);
end
if ischar(Func2)
    Func2 = str2func(Func2);
end

Out1 = cell(1, Nout);
Out2 = cell(1, Nout);

tic
for k = 1:Iter
    [Out1{:}] = Func1(varargin{:});
end
Time1 = toc;

tic
for k = 1:Iter
    [Out2{:}] = Func2(varargin{:});
end
Time2 = toc;

AssertTest = zeros(1, Nout, 'logical');
for j = 1:Nout
    AssertTest(j) = isequal(Out1{j}, Out2{j});
end

Summary.Time1 = Time1;
Summary.Time2 = Time2;
Summary.Ratio = Time1/Time2;
Summary.AssertTest = AssertTest;