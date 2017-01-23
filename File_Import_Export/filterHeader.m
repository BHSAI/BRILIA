%filterHeader will check the cell array for presence of a header. A header
%is simply a row with nothing but character values. If detected, this will
%trim off the header, returning the data portion and header portion
%separately. 
%
%  [Data, Header] = filterHeader(CellData) will remove Header rows from the
%  CellData, and return it as Data. Header will also be saved.
%
%  [Data, Header] = filterHeader(CellData,N) will remove N rows from the
%  CellData, and return it as Data. Header will also be saved.
function varargout = filterHeader(CellData,varargin)

if isempty(varargin)
    IsHeader = 1; %First assume there is a header
    for j = 1:size(CellData,1)
        %Check to see if all entries in rows are all characters
        for h = 1:size(CellData,2)
            if isnumeric(CellData{j,h}) && (isnan(CellData{j,h}) == 0)
                IsHeader = 0;
                break;
            end
        end
        if IsHeader == 0
            break
        end
    end
    %Possible that you get a string array. Need to then assume 1st row is it.
    if j == size(CellData,1)
        j = 1;
    else
        j = j-1; %The header ended 1 row before.
    end
else
    j = varargin{1};
end

%Return the Data and Header separately.
if nargout >= 1
    varargout{1} = CellData(j+1:end,:);
    if nargout == 2
        Header = CellData(1:j,:);
        for j = 1:length(Header)
            if  sum(isnan(Header{j}))>0 || sum(isempty(Header{j}))>0; 
                Header{j} = '';
            else
                Header{j} = strrep(Header{j},'"','');
            end
        end
        varargout{2} = Header; %CellData(1:j,:);
    end
end
