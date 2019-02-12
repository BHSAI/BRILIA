%applyall will apply a function for multiple data columns and returns a
%group-treated column based on the Grouping vector. This is used to take
%mean, std, median, etc of data set by grouping columns together.
%
%  INPUT
%    Func ['mean', 'mode', 'median', 'std', 'min', 'max', 'var']: function 
%      string or handle to use for doing column-wise calculations
%    Data: MxN data set where N columns are for length(unique(N)) groups
%    Grouping: 1xN vector assign each column of Data to a treatement group
%      number 1, 2, 3, etc. 
%
%  OUTPUT
%    Out: Mxlength(unique(Grouping)) vector of results
%
%  EXAMPLE
%    Data = [rand(3, 5) 10*rand(3, 4) -10*rand(3,6)];   %5-, 4-, 6-col data 
%    Grouping = repelem([1 2 3], [5 4 6]); %[1 1 1 1 1 2 2 2 2 3 3 3 3 3 3]
%    Mean = applyall('mean', Data, Grouping); %mean of 1's, 2's, and 3's
%
function [Out, UnqGrp] = applyall(Func, Data, Grouping)
if ischar(Func)
    Func = str2func(Func);
elseif ~ishandle(Func)
    error('%s: Func must be a string or function handle.', mfilename);
end

UnqGrp = unique(Grouping);
Out = zeros(size(Data, 1), length(UnqGrp));
StrFunc = func2str(Func);

if contains(StrFunc, {'mean', 'median', 'mode'})
    for g = 1:length(UnqGrp)
        Out(:, g) = Func(Data(:, UnqGrp(g) == Grouping), 2);
    end
elseif contains(StrFunc, {'std', 'min', 'max', 'var'})
    for g = 1:length(UnqGrp)
        Out(:, g) = Func(Data(:, UnqGrp(g) == Grouping), [], 2);
    end
else
    for g = 1:length(UnqGrp)
        Loc = UnqGrp(g) == Grouping;
        for k = 1:size(Data, 1)
            Out(k, g) = Func(Data(k, Loc));
        end
    end
end