%findCell will look in cell of A for any match to cell of B.
%
%  Aidx = findCell(ListA,ListB) will search ListA that matches with
%  ListB, where B is either a character or cell of string, number, NaN, or
%  []. If there is no match, returns 0.
%
%  Aidx = findCell(ListA,ListB,'anycase') will do a cell string
%  search, but for any case match.
%
%  EX for string search: 
%    ListA = {'a' [] 'C' 'c' 'D' '' NaN 3 4}
%
%    ListB = {'c' 'd'}
%    Aidx = findCell(ListA,ListB)
%        Aidx = 4
%
%  EX for string any case search: 
%    ListB = {'c' 'd'}
%    Aidx = findCell(ListA,ListB,'anycase')
%        Aidx = 3
%               4
%               5
%
%  EX for number search: 
%    ListB = {3 4}
%    Aidx = findCell(ListA,ListB)
%        Aidx = 7
%               8
%
%  EX for NaN cell search: 
%    ListB = NaN
%    Aidx = findCell(ListA,ListB)
%        Aidx = 7
%
%  EX for empty cell search: 
%    ListB = [] 
%    Aidx = findCell(ListA,ListB)
%        Aidx = 2
%               6

function Aidx = findCell(ListA,ListB,varargin)
%Make sure ListA and ListB are all cells, and a list.
if ~iscell(ListA)
    ListA = {ListA};
end
if ~iscell(ListB)
    ListB = {ListB};
end
if min(size(ListA)) > 1
    error('ListA must be a 1xN or Nx1 cell matrix');
end
if min(size(ListB)) > 1
    error('ListB must be a 1xN or Nx1 cell matrix');
end

%if you have an option for 'anycase', then assume both inputs are str, and
%set to lower cases.
if ~isempty(varargin)
    if strcmpi(varargin{1},'anycase')
        for a = 1:length(ListA)
            if ischar(ListA{a})
                ListA{a} = lower(ListA{a});
            end
        end
        for b = 1:length(ListB)
            if ischar(ListB{b})
                ListB{b} = lower(ListB{b});
            end
        end
    end
end

%Determine the class of B
if ischar(ListB{1})
    ClassType = 'char';
elseif isnumeric(ListB{1});
    ClassType = 'numeric';
    if isnan(ListB{1})
        ClassType = 'nan';
    elseif isempty(ListB{1});
        ClassType = 'empty';
    end
else
    disp('findCell error: unknown class type of ListA. Defaulting to char class')
    ClassType = 'char';
end

%Perform search appropriate for each
switch ClassType
    case 'char'
        for j = 1:length(ListA)
            if ~ischar(ListA{j})
                ListA{j} = '';
            end
        end
        Aidx = find(ismember(ListA,ListB))';
    case 'numeric'
        for j = 1:length(ListA)
            if ~isnumeric(ListA{j})
                ListA{j} = NaN;
            end
        end
        [~,Aidx,~] = intersect(cell2mat(ListA),cell2mat(ListB));
    case 'nan'
        Aidx = zeros(length(ListA),1,'logical');
        for j = 1:length(ListA)
            if isnan(ListA{j})
                Aidx(j) = 1;
            end
        end
        Aidx = find(Aidx);
    case 'empty'
        Aidx = zeros(length(ListA),1,'logical');
        for j = 1:length(ListA)
            if isempty(ListA{j})
                Aidx(j) = 1;
            end
        end
        Aidx = find(Aidx);
end

%Ensure Aidx is 0 for empty.
if isempty(Aidx)
    Aidx = 0;
end