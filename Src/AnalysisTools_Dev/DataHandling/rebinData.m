%rebinData will redo the data binning fo a M x (1+N) cell array, where the
%1st column is the string/numeric fields and the rest are the frequencies.
%
%  R = rebinData(FreqCell)
%
%  R = rebinData(FreqCell, Edges)
%
%  INPUT
%    FreqCell: M x (N+1) cell of frequencies, where FreqCell(:, 1) are the unique, sorted
%      labels from all distributions, which are the outputs from countData
%      or countData.
%    Edges: vector for rebinning numeric fields of FreqCell. 
%
%  OUTPUT
%    R: rebinned data of FreqCell based
%
%  NOTE
%    Edges MUST be specified in order to rebin data with numeric fields.
%
%    If FreqCell(:, 1) are string, rebinData will get the unique char cells and
%    add up the frequencies. This is useful if you want to regroup the x
%    labels.
%
%  EXAMPLE
%    FreqCell = {'a'  1    2  3;
%            'b'  1    2  3;
%            'c'  1    2  3;
%            3    3.1  4  5;
%            3.2  3.2  2  1;
%            4    1.3  3  1;
%            4.2  1.3  2  3};
%    Edges = [0:5];
%    FreqCell(1:2, 1) = {'a'};
%    R = rebinData(FreqCell, Edges);
%    R =
%         'a'    [     2]    [4]    [6]
%         'c'    [     1]    [2]    [3]
%         [0]    [     0]    [0]    [0]
%         [1]    [     0]    [0]    [0]
%         [2]    [     0]    [0]    [0]
%         [3]    [6.3000]    [6]    [6]
%         [4]    [2.6000]    [5]    [4]
%         [5]    [     0]    [0]    [0]
%
function ReFreqCell = rebinData(FreqCell, varargin)
if ~iscell(FreqCell)
    FreqCell = num2cell(FreqCell);
end
StrLoc = cellfun('isclass', FreqCell(:, 1), 'char'); 
StrData = recountStrData(FreqCell( StrLoc, :));
NumData = recountNumData(FreqCell(~StrLoc, :), varargin{:});
ReFreqCell = vertcat(StrData, NumData);

%Recount string field frequencies, returning only their bin and frequency
function StrData = recountStrData(Data)
Name = Data(:, 1);
Freq = cell2mat(Data(:, 2:end));
[ReName, ~, ~, ReNameIdx] = unique2(Name);
ReFreq = zeros(numel(ReName), size(Freq, 2));
for j = 1:numel(ReName)
    ReFreq(j, :) = sum(Freq(ReNameIdx{j}, :), 1);
end
StrData = [ReName num2cell(ReFreq)];

%Recount numeric field frequencies, returning only their bin and frequency
function NumData = recountNumData(Data, varargin)
if isempty(varargin)
    NumData = Data;
    return
end
Name = cell2mat(Data(:, 1));
Freq = cell2mat(Data(:, 2:end));
[BinNum, ReName] = discretize(Name, varargin{:});
[UnqBinNum, ~, ~, UnqBinIdx] = unique2(BinNum);
ReFreq = zeros(numel(ReName), size(Freq, 2));
for j = 1:numel(UnqBinNum)
    if isnan(UnqBinNum(j)); continue; end
    ReFreq(UnqBinNum(j), :) = sum(Freq(UnqBinIdx{j}, :), 1);
end
NumData = num2cell([ReName(:) ReFreq]);