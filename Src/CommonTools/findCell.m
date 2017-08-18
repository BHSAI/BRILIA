%findCell will look in ListA for any match ListB, returning the index
%of any matches in ListA. ListA should be the main cell matrix your are
%looking for the index value.
%
%  Aidx = findCell(SearchList,QueryListB)
%
%  Aidx = findCell(SearchQueryList, QueryListB, Param, Value, ...) will allow certain
%
%  INPUT
%    SearchQueryList: 1xN or Nx1 QueryList of string or nubmers
%    QueryListB: 1xN or Nx1 QueryList of contents you are searching for in SearchQueryList
%    Param, Value:
%       PARAM       VALUE               DESCRIPTION
%       ----------  ------------------  -------------------------------
%       MatchCase   'exact' or 'any'    Exact case or any case match    
%       MatchWord   'all' or 'partial'  full word or partial word ma
%
%  OUTPUT
%    Aidx: The location in SearchQueryList that matches with any item in QueryListB.
%      Returns 0 for no match. Results are in the order of QueryListB.
%
%  EXAMPLE 
%      SearchQueryList = {'a' [] 'C' 'c'  'D' '' NaN 3 4 'cat' 'dog'}
%
%    Case1) For string searches
%      QueryListB = {'c' 'd'}
%      Aidx = findCell(SearchQueryList, QueryListB)
%          Aidx = 4
%
%    Case2) For string any case searches
%      QueryListB = {'c' 'd'}
%      Aidx = findCell(SearchQueryList, QueryListB, 'MatchCase', 'any', 'MatchWord', 'all')
%          Aidx = 3
%                 4
%                 5
%
%    Case3) For string partial match searches
%      QueryListB = {'c' 'd'}
%      Aidx = findCell(SearchQueryList, QueryListB, 'MatchCase', 'exact', 'MatchWord', 'partial')
%          Aidx = 4
%                 10
%                 11
%
%    Case4) For number searches, and QueryListB order is different than SearchQueryList
%    QueryListB = {4 3}
%    Aidx = findCell(SearchQueryList, QueryListB)
%          Aidx = 8
%                 7
%
%    Case5) For NaN cell search: 
%      QueryListB = NaN
%      Aidx = findCell(SearchQueryList, QueryListB)
%          Aidx = 7
%
%    Case6) For empty cell search: 
%      QueryListB = [] 
%      Aidx = findCell(SearchQueryList, QueryListB)
%          Aidx = 2
%                 6

function Aidx = findCell(SearchQueryList, QueryListB, varargin)
%Parse parameter inputs
P = inputParser;
addParameter(P, 'MatchCase', 'exact', @ischar);
addParameter(P, 'MatchWord', 'all', @ischar);
parse(P, varargin{:});

%Make sure SearchQueryList and QueryListB are all cells, and a QueryList.
if ~iscell(SearchQueryList)
    SearchQueryList = {SearchQueryList};
end
if ~iscell(QueryListB)
    QueryListB = {QueryListB};
end
if min(size(SearchQueryList)) > 1
    error('SearchQueryList must be a 1xN or Nx1 cell matrix');
end
if min(size(QueryListB)) > 1
    error('QueryListB must be a 1xN or Nx1 cell matrix');
end

%If you have an option for 'anycase', then assume both inputs are str, and
%set to lower cases.
if strcmpi(P.Results.MatchCase, 'Any')
    for a = 1:length(SearchQueryList)
        if ischar(SearchQueryList{a})
            SearchQueryList{a} = lower(SearchQueryList{a});
        end
    end
    for b = 1:length(QueryListB)
        if ischar(QueryListB{b})
            QueryListB{b} = lower(QueryListB{b});
        end
    end
end

%Determine the class of B
if ischar(QueryListB{1})
    ClassType = 'char';
elseif isnumeric(QueryListB{1});
    ClassType = 'numeric';
    if isnan(QueryListB{1})
        ClassType = 'nan';
    elseif isempty(QueryListB{1});
        ClassType = 'empty';
    end
else
    disp('findCell error: unknown class type of SearchQueryList. Defaulting to char class')
    ClassType = 'char';
end

%Perform search appropriate for each
switch ClassType
    case 'char'
        for j = 1:length(SearchQueryList)
            if ~ischar(SearchQueryList{j})
                SearchQueryList{j} = '';
            end
        end
        if strcmpi(P.Results.MatchWord, 'All')
            [Aidxb, Bidx] = ismember(SearchQueryList, QueryListB);
            [~, SortIdx] = sort(Bidx(Aidxb));
            Aidx = find(Aidxb)';
            Aidx = Aidx(SortIdx);
        else
            Bidx = zeros(length(SearchQueryList), 1, 'logical');
            for b = 1:length(QueryListB)
                for k = 1:length(SearchQueryList)
                    if ~isempty(regexp(SearchQueryList{k}, QueryListB{b}, 'once'))
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
        for j = 1:length(SearchQueryList)
            if ~isnumeric(SearchQueryList{j})
                SearchQueryList{j} = NaN;
            end
        end
        [~, Aidx, Bidx] = intersect(cell2mat(SearchQueryList), cell2mat(QueryListB));
        [~, SortIdx] = sort(Bidx);
        Aidx = Aidx(SortIdx);
    case 'nan'
        Aidx = zeros(length(SearchQueryList), 1, 'logical');
        for j = 1:length(SearchQueryList)
            if isnan(SearchQueryList{j})
                Aidx(j) = 1;
            end
        end
        Aidx = find(Aidx);
    case 'empty'
        Aidx = zeros(length(SearchQueryList), 1, 'logical');
        for j = 1:length(SearchQueryList)
            if isempty(SearchQueryList{j})
                Aidx(j) = 1;
            end
        end
        Aidx = find(Aidx);
end

%Ensure Aidx is 0 for empty.
if isempty(Aidx)
    Aidx = 0;
end
