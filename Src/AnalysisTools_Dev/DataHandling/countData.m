%countData will take in a cell array of string and/or numbers and return
%the frequency of each cell in a Mx2 cell output. If you have have a "|" or
%a "," divider, then it will split up the cell and assign a weight to it.
%
%  FreqCell = countData(C)
%
%  FreqCell = countData(C, Weight)
%
%  INPUT
%    C: Mx1 cell array
%    Weight: Mx1 matrix of weight per each C
%
%  OUTPUT
%    FreqCell: Ux2 cell where 
%      Col 1 - string or number
%      Col 2 - frequency
%
%  EXAMPLE
%    C = {'a|b' 'b' 'c' 'c' 'c' 1 2 2 2 3 3 4}
%    FreqCell = countData(C)
%    FreqCell =
%         [1]    [     1]
%         [2]    [     3]
%         [3]    [     2]
%         [4]    [     1]
%         'a'    [0.5000]
%         'b'    [1.5000]
%         'c'    [     3]
%
%    C = {'2IGHD2-3*01';
%         '3IGHD5-5*01';
%         '2IGHD5-2*01|2IGHD5-3*01|2IGHD5-4*01|2IGHD5-6*01';
%         '3IGHD4-1*01';
%         '1rIGHD2-5*01|1rIGHD2-6*01';
%         '2IGHD2-4*01';
%         '1IGHD4-1*01';
%         '3IGHD2-3*01';
%         '1IGHD2-4*01';
%         '2rIGHD3-1*01'};
%    FreqCell = countData(C)
%    FreqCell =
%         '1IGHD2-4*01'     [     1]
%         '1IGHD4-1*01'     [     1]
%         '1rIGHD2-5*01'    [0.5000]
%         '1rIGHD2-6*01'    [0.5000]
%         '2IGHD2-3*01'     [     1]
%         '2IGHD2-4*01'     [     1]
%         '2IGHD5-2*01'     [0.2500]
%         '2IGHD5-3*01'     [0.2500]
%         '2IGHD5-4*01'     [0.2500]
%         '2IGHD5-6*01'     [0.2500]
%         '2rIGHD3-1*01'    [     1]
%         '3IGHD2-3*01'     [     1]
%         '3IGHD4-1*01'     [     1]
%         '3IGHD5-5*01'     [     1]
%
%  NOTE
%    If C contains a string with multiple options (Ex: C{c} = 'a|b|c|'),
%    then the Weight of that will be split by the number of elements.
%
function FreqCell = countData(C, Weight)
if nargin == 1
    Weight = ones(size(C));
elseif nargin == 2 && iscell(Weight)
    Weight(cellfun('isempty', Weight)) = {0};
end
if numel(C) ~= numel(Weight)
    error('%s: Number of cell must be same as number of weights', mfilename);
end
if ~iscell(C)
    C = num2cell(C);
end
if ~iscell(Weight)
    Weight = num2cell(Weight);
end

%Split up a dividing type.
NumElem = 0;
for k = 1:numel(C)
    if ischar(C{k}) && any(contains(C{k}, {'|', ','}))
        C{k} = strsplit(C{k}, {',', '|'})';
        Weight{k} = repelem(Weight{k}/length(C{k}), length(C{k}), 1);
        NumElem = NumElem + length(C{k});
    else
        C{k} = C(k);
        NumElem = NumElem + 1;
    end
end

C = vertcat(C{:});
Weight = vertcat(Weight{:});
WIdx = 1:length(Weight);

%Do the number data first
NumLoc = cellfun(@isnumeric, C);
WNumIdx = WIdx(NumLoc);
[UnqNum, ~, IdxNum] = unique(cell2mat(C(NumLoc)), 'stable');
FreqNum = cell(length(UnqNum), 2);
for k = 1:length(UnqNum)
    FreqNum(k, :) = {UnqNum(k) sum(Weight(WNumIdx(IdxNum == k)))};
end

%Do the text data next
WTxtIdx = WIdx(~NumLoc);
[UnqTxt, ~, IdxTxt] = unique(C(~NumLoc));
FreqTxt = cell(length(UnqTxt), 2);
for k = 1:length(UnqTxt)
    FreqTxt(k, :) = {UnqTxt{k} sum(Weight(WTxtIdx(IdxTxt == k)))};
end

FreqCell = vertcat(FreqTxt, FreqNum);
[~, SortIdx] = sort2(FreqCell(:, 1));
FreqCell = FreqCell(SortIdx, :);