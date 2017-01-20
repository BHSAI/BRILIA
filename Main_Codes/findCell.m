%findCell will look in ListA for any match ListB, returning the index
%of any matches in ListA.
%
%  Aidx = findCell(ListA,ListB)
%
%  Aidx = findCell(ListA,ListB,Param,Value,...) will allow certain
%
%  INPUT
%    ListA: 1xN or Nx1 list of string or nubmers
%    ListB: 1xN or Nx1 list of contents you are searching for in ListA
%    Param,Value:
%       PARAM       VALUE               DESCRIPTION
%       ----------  ------------------  -------------------------------
%       MatchCase   'exact' or 'any'    Exact case or any case match    
%       MatchWord   'all' or 'partial'  full word or partial word ma
%
%  OUTPUT
%    Aidx: The location in ListA that matches with any item in ListB.
%      Returns 0 for no match.
%
%  EXAMPLE 
%      ListA = {'a' [] 'C' 'c'  'D' '' NaN 3 4 'cat' 'dog'}
%
%    Case1) For string searches
%      ListB = {'c' 'd'}
%      Aidx = findCell(ListA,ListB)
%          Aidx = 4
%
%    Case2) For string any case searches
%      ListB = {'c' 'd'}
%      Aidx = findCell(ListA,ListB,'MatchCase','any','MatchWord','all')
%          Aidx = 3
%                 4
%                 5
%
%    Case3) For string partial match searches
%      ListB = {'c' 'd'}
%      Aidx = findCell(ListA,ListB,'MatchCase','exact','MatchWord','partial')
%          Aidx = 4
%                 10
%                 11
%
%    Case4) For number searches 
%    ListB = {3 4}
%    Aidx = findCell(ListA,ListB)
%          Aidx = 7
%                 8
%
%    Case5) For NaN cell search: 
%      ListB = NaN
%      Aidx = findCell(ListA,ListB)
%          Aidx = 7
%
%    Case6) For empty cell search: 
%      ListB = [] 
%      Aidx = findCell(ListA,ListB)
%          Aidx = 2
%                 6

function Aidx = findCell(ListA,ListB,varargin)
%Parse parameter inputs
P = inputParser;
addParameter(P,'MatchCase','exact',@ischar);
addParameter(P,'MatchWord','all',@ischar);
parse(P,varargin{:});

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

%If you have an option for 'anycase', then assume both inputs are str, and
%set to lower cases.
if strcmpi(P.Results.MatchCase,'Any')
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
        if strcmpi(P.Results.MatchWord,'All')
            Aidx = find(ismember(ListA,ListB))';
        else
            Aidx = zeros(length(ListA),1,'logical');
            Bphrase = ListB{1};
            for b = 2:length(ListB)
                Bphrase = [Bphrase '|' ListB{b}];
            end
            for k = 1:length(ListA)
                if ~isempty(regexp(ListA{k},Bphrase,'once'))
                    Aidx(k) = 1;
                end
            end
            Aidx = find(Aidx);
        end
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
