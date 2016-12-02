%findHeader will look at the header cells and look for the location of a
%certain header query. Finds the first one only default
%
%  HeaderLoc = findHeader(Header,Query) this is default lookup. WILL NOT
%  LOOK if there are same hits. Returns a size(Query) matrix of locations
%  hit in Header index.
%
%  EX: 
%    Header = {'a' 'b' 'c' 'd'}
%    HeaderLoc = findHeader(Header,{'c' 'd'})
%      HeaderLoc = [3 4]

function HeaderLoc = findHeader(HeaderRef,HeaderCheck,varargin)
%Make sure HeaderCheck is a cell
if ~iscell(HeaderCheck) && ischar(HeaderCheck)
    HeaderCheck = {HeaderCheck};
end

CheckAll = 'first';
if ~isempty(varargin)
    CheckAll = varargin{1};
    if ~strcmpi(CheckAll,'all') && ~strcmpi(CheckAll,'first')
        error('3rd option should be ''all'' or ''first''');
    end        
end

%Look for the Header Loc.
if strcmpi(CheckAll,'first')
    HeaderLoc = zeros(1,length(HeaderCheck));
    for k = 1:length(HeaderCheck);
        for j = 1:length(HeaderRef)
            if strcmpi(HeaderRef{j},HeaderCheck{k})
                HeaderLoc(k) = j;
                break
            end
        end
    end
else
    if length(HeaderCheck) > 1
        error('Cannot have multiple HeaderCheck when using ''all'' option');
    end
    HeaderLoc = zeros(1,length(HeaderRef));
    for j = 1:length(HeaderRef)
        if strcmpi(HeaderRef{j},HeaderCheck{1})
            HeaderLoc(j) = j;
        end
    end
    HeaderLoc(HeaderLoc == 0) = [];
end