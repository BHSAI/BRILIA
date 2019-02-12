%testall is a wrapper for doing multiple hypothesis testing on sets of data
%beloning to N groups, in which N*(N-1)/2 pairwise comparisons will be
%performed between the data groups.
%
%  [Out1, Out2, Out3, Out4, Comb] = testall(Func, Data, Grouping, varargin)
%
%  INPUT
%    Func: function string or handle to use for doing the test
%    Data: MxN data set for M parameters and N1 subjects
%    Grouping: 1xN vector assign each column of Data to a treatement group
%      number 1, 2, 3, etc. 
%    varargin: any additional parameters required for the Func used
%
%  OUTPUT
%   OutN: output 1-4 that is normally returned from the function used.
%   Comb: the 5th output is simply nchoosek(unique(Grouping), 2)' that
%   shows the ith and jth group that are compared against using the
%   Func(Grp_i, Grp_j) command.
%
%  EXAMPLE
%    Data = [rand(3, 5) 10*rand(3, 5) -10*rand(3,6) rand(3,6)]
%    Grouping = repelem([1 2 3 4], [5 5 6 6]);
%    [H, P, CI, Stat, Comb] = testall('ttest2', Data, Grouping, 'alpha', 0.05)
%    H =
%          1     1     0     1     1     1
%          1     1     0     1     1     1
%          1     1     0     1     1     1
% 
%    P =
%         0.0416    0.0120    0.3524    0.0026    0.0207    0.0075
%         0.0003    0.0215    0.0808    0.0005    0.0001    0.0089
%         0.0093    0.0016    0.0738    0.0003    0.0026    0.0011
% 
%    CI =
%       3×6 cell array of 1x2 double
% 
%    Stat =
%       3x6 cell array of 1x1 struct
% 
%    Comb =
%          1     1     1     2     2     3
%          2     3     4     3     4     4
%
%    [H, P, CI, DXY, Comb] = testall('donDiffTest', Data, Grouping, 'alpha', 0.05)
%    H =
%          0     1     0     1     0     1
%          1     1     0     1     1     1
%          0     1     0     1     1     1
% 
%    P =
%         0.4000    0.0000    0.6000    0.0000    0.3333    0.0000
%         0.0000    0.0000    0.3333    0.0000    0.0000    0.0000
%         0.2400    0.0000    0.4667    0.0000    0.0000    0.0000
% 
%    CI =
%       3×6 cell array of 1x2 double
% 
%    DXY =
%       3x6 cell array of 1xN double
% 
%    Comb =
%          1     1     1     2     2     3
%          2     3     4     3     4     4
     
function varargout = testall(Func, Data, Grouping, varargin)
if ischar(Func)
    Func = str2func(Func);
elseif ~ishandl(Func)
    error('%s: Func must be a string or function handle.', mfilename);
end

[UnqGrp, ~, UnqIdx] = unique(Grouping, 'stable');
UnqCellIdx = cell(1, length(UnqGrp));
for y = 1:length(UnqGrp)
    UnqCellIdx{y} = find(UnqIdx == y);
end

Nout = nargout;
varargout = cell(1, Nout);
CombGrp = nchoosek(1:length(UnqGrp), 2)';
AllOut = cell(size(Data, 1), size(CombGrp, 1));
if Nout > 4
    Nout = 4;
    varargout{5} = CombGrp;
end

for k = 1:size(Data, 1)
    for c = 1:size(CombGrp, 2)
        [AllOut{k, c}{1:Nout}] = Func(Data(k, UnqCellIdx{CombGrp(1, c)}), Data(k, UnqCellIdx{CombGrp(2, c)}), varargin{:});
    end
end

for n = 1:Nout
    Temp = cellfun(@(x) x{n}, AllOut, 'unif', false);
    if all(all(cellfun(@(x) isnumeric(x) && numel(x)==1, Temp)))
        Temp = cell2mat(Temp);
    end
    varargout{n} = Temp;
end