%sort2 will sort a cell array of text and number according to
%grouped chars and num. 
%
%  [SortedList, SortIdx] = sortTextAnNum(List, varargin)
%
%  INPUT
%    List: cell array of string and num
%    varargin: inputs that can be used by sort.m
%
%  OUTPUT
%    SortedList: sorted list
%    SortIdx: sort index order so S = List(Idx)
%
%  EXAMPLE
%    List = {'ABC-1'; 'ABB-1'; 'ABB-12'; 'ABB-4'};
%    [SortedList, SortIdx] = sort2(List)
%    SortedList =
%         'ABB-1'
%         'ABB-4'
%         'ABB-12'
%         'ABC-1'
% 
%     SortIdx =
%          2
%          4
%          3
%          1
%
%    List = {'1-ABC'; '1-ABB-1'; '1-ABB-11';'12-ABB'; '4-ABB'};
%    [SortedList, SortIdx] = sort2(List)
% 
%    SortedList =
%         '1-ABB-1'
%         '1-ABB-11'
%         '1-ABC'
%         '4-ABB'
%         '12-ABB'
%     SortIdx =
%          2
%          3
%          1
%          5
%          4
function [SortedList, SortIdx] = sort2(List, varargin)
if isnumeric(List) %Sorting only numbers
    [SortedList, SortIdx] = sort(List, varargin{:});
    return
elseif ischar(List) %Sorting a single char line (do nothing)
    SortedList = List;
    SortIdx = 1;
    return
elseif ~iscell(List)
    error('%s: 1st input should be a matrix or cell of string and/or num.', mfilename);
end

%From here one, dealing with sorting cell arrays
if all(cellfun('isclass', List, 'double')) %Soring a cell of only numbers
    [SortedList, SortIdx] = sort(cell2mat(List), varargin{:});
    SortedList = num2cell(SortedList);
    return
end

%Sorting char and/or numbers in a cell of char and/or numbers
SplitCell = cell(numel(List), 1);
for j = 1:numel(List)
    SplitCell{j} = splitTextNum(List{j});
end
MaxNum = max(cellfun('length', SplitCell));

TempList = num2cell(zeros(numel(List), MaxNum));
for j = 1:numel(List)
    TempList(j, 1:length(SplitCell{j})) = SplitCell{j};
end

NumLoc = cellfun('isclass', TempList, 'double');
LeadZeros = max(ceil(log10(cell2mat(TempList(NumLoc)))));
if LeadZeros <= 0
    NumForm = '%f';
else
    NumForm = sprintf('%%0%dd', LeadZeros);
end
TempList(NumLoc) = cellfun(@(x) sprintf(NumForm, round(x)), TempList(NumLoc), 'unif', false);
SortList = cell(numel(List), 1);
for k = 1:length(SortList)
    SortList{k} = [TempList{k, :}];
end
if ~isempty(varargin)
    warning('%s: DIM and MODE not supported for cell of mixed str and num content.', mfilename);
end
[~, SortIdx] = sort(SortList);
SortedList = List(SortIdx);

function S = splitTextNum(TextNum)
if isnumeric(TextNum)
    S = {TextNum};
    return
end
NumLoc = isstrprop(TextNum, 'digit');
DecimalIdx = regexp(TextNum, '\d\.\d')+1;
NumLoc(DecimalIdx) = 1;
RegLoc = labelRegionMEX(NumLoc);
S = cell(1, max(RegLoc));
for k = 1:max(RegLoc)
    Loc = RegLoc == k;
    if any(NumLoc(Loc))
        try
            S{k} = convStr2NumMEX(TextNum(Loc));
        catch
            S{k} = TextNum(Loc);
        end
    else
        S{k} = TextNum(Loc);
    end
end
