%findCell will look in ListA for any match ListB, returning the index
%of any matches in ListA. ListA should be the main cell matrix your are
%looking for the index value.
%
%  Aidx = findCell(SearchList, QueryList)
%
%  Aidx = findCell(SearchList, QueryList, Param, Value, ...) will allow certain
%
%  INPUT
%    SearchList: 1xN or Nx1 QueryList of string or nubmers
%    QueryList: 1xN or Nx1 QueryList of contents you are searching for in SearchList
%    Param, Value:
%       PARAM       VALUE               DESCRIPTION
%       ----------  ------------------  -------------------------------
%       MatchCase   'exact' or 'any'    Exact case or any case match    
%       MatchWord   'all' or 'partial'  full word or partial word ma
%
%  OUTPUT
%    Aidx: The location in SearchList that matches with any item in QueryList.
%      Returns 0 for no match. Results are in the order of QueryList.
%
%  EXAMPLE 
%      SearchList = {'a' [] 'C' 'c'  'D' '' NaN 3 4 'cat' 'dog'}
%
%    Case1) For string searches
%      QueryList = {'c' 'd'}
%      Aidx = findCell(SearchList, QueryList)
%          Aidx = 4
%
%    Case2) For string any case searches
%      QueryList = {'c' 'd'}
%      Aidx = findCell(SearchList, QueryList, 'MatchCase', 'any', 'MatchWord', 'all')
%          Aidx = 3
%                 4
%                 5
%
%    Case3) For string partial match searches
%      QueryList = {'c' 'd'}
%      Aidx = findCell(SearchList, QueryList, 'MatchCase', 'exact', 'MatchWord', 'partial')
%          Aidx = 4
%                 10
%                 11
%
%    Case4) For number searches, and QueryList order is different than SearchList
%    QueryList = {4 3}
%    Aidx = findCell(SearchList, QueryList)
%          Aidx = 8
%                 7
%
%    Case5) For NaN cell search: 
%      QueryList = NaN
%      Aidx = findCell(SearchList, QueryList)
%          Aidx = 7
%
%    Case6) For empty cell search: 
%      QueryList = [] 
%      Aidx = findCell(SearchList, QueryList)
%          Aidx = 2
%                 6

function Aidx = findCell(SearchList, QueryList, varargin)
P = inputParser;
addParameter(P, 'MatchCase', 'exact', @ischar);
addParameter(P, 'MatchWord', 'all', @ischar);
parse(P, varargin{:});

%Make sure SearchList and QueryList are all cells, and a QueryList.
if ~iscell(SearchList)
    SearchList = {SearchList};
end
if ~iscell(QueryList)
    QueryList = {QueryList};
end
if min(size(SearchList)) > 1
    error('SearchList must be a 1xN or Nx1 cell matrix');
end
if min(size(QueryList)) > 1
    error('QueryList must be a 1xN or Nx1 cell matrix');
end

%If you have an option for 'anycase', then assume both inputs are str, and
%set to lower cases.
if strcmpi(P.Results.MatchCase, 'Any')
    for a = 1:length(SearchList)
        if ischar(SearchList{a})
            SearchList{a} = lower(SearchList{a});
        end
    end
    for b = 1:length(QueryList)
        if ischar(QueryList{b})
            QueryList{b} = lower(QueryList{b});
        end
    end
end

%Determine the class of B
if ischar(QueryList{1})
    ClassType = 'char';
elseif isnumeric(QueryList{1})
    ClassType = 'numeric';
    if isnan(QueryList{1})
        ClassType = 'nan';
    elseif isempty(QueryList{1})
        ClassType = 'empty';
    end
else
    warning('%s: unknown class type of SearchList. Defaulting to char class.', mfilename)
    ClassType = 'char';
end

%Perform search appropriate for each
switch ClassType
    case 'char'
        for j = 1:length(SearchList)
            if ~ischar(SearchList{j})
                SearchList{j} = '';
            end
        end
        if strcmpi(P.Results.MatchWord, 'All')
            [Aidxb, Bidx] = ismember(SearchList, QueryList);
            [~, SortIdx] = sort(Bidx(Aidxb));
            Aidx = find(Aidxb)';
            Aidx = Aidx(SortIdx);
        else
            Bidx = zeros(length(SearchList), 1, 'logical');
            for b = 1:length(QueryList)
                for k = 1:length(SearchList)
                    if ~isempty(regexp(SearchList{k}, QueryList{b}, 'once'))
                        Bidx(k) = b;
                    end
                end
            end            
            Aidxb = Bidx > 0;
            [~, SortIdx] = sort(Bidx(Aidxb));
            Aidx = find(Aidxb)';
            Aidx = Aidx(SortIdx);
        end
    case 'numeric'
        for j = 1:length(SearchList)
            if ~isnumeric(SearchList{j})
                SearchList{j} = NaN;
            end
        end
        [~, Aidx, Bidx] = intersect(cell2mat(SearchList), cell2mat(QueryList));
        [~, SortIdx] = sort(Bidx);
        Aidx = Aidx(SortIdx);
    case 'nan'
        Aidx = zeros(length(SearchList), 1, 'logical');
        for j = 1:length(SearchList)
            if isnan(SearchList{j})
                Aidx(j) = 1;
            end
        end
        Aidx = find(Aidx);
    case 'empty'
        Aidx = zeros(length(SearchList), 1, 'logical');
        for j = 1:length(SearchList)
            if isempty(SearchList{j})
                Aidx(j) = 1;
            end
        end
        Aidx = find(Aidx);
end

%Ensure Aidx is 0 for empty.
if isempty(Aidx)
    Aidx = 0;
end
