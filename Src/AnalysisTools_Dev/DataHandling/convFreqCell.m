%convFreqTable will convert a Mx(N+1) cell into a Mx1 header cell and MxN
%matrix cell, and vice versa.
%
%  [Category, Freq] = convFreqCell(FreqCell)
%
%  FreqCell = convFreqCell(Category, Freq)
%
%  INPUT/OUTPUT
%    FreqCell: Mx(N+1) cell where 1st column is the cateogry and the other
%       N columns are frequency for N distributions or subjects.
%    Category: Mx1 cell of categories, 1st column of FreqCell
%    Freq: MxN frequency matrix
%
%  EXAMPLE
%    FreqCell = {'a' 20 30; 
%                 'b' 30 50}
%    [Category, Freq] = convFreqCell(FreqCell)
%    Category =
%            'a'
%            'b'
%    Freq =
%            20    30
%            30    50
%
%    FreqCell = convFreqCell(Category, Freq)
%    FreqCell =
%         'a'    [20]    [30]
%         'b'    [30]    [50]
function varargout = convFreqCell(varargin)
if nargin == 1 && iscell(varargin{1})
    if nargout ~= 2
        error('%s: Must have 2 output when using 1 input.', mfilename);
    end
    varargout{1} = varargin{1}(:, 1);
    varargout{2} = cell2mat(varargin{1}(:, 2:end));
elseif nargin == 2 && iscell(varargin{1}) && isnumeric(varargin{2})
    if numel(varargin{1}) ~= size(varargin{2}, 1)
        error('%s: Number of Input1 elements (%d) not same as Input2 rows (%d). ', mfilename, numel(varargin{1}), size(varargin{2}, 1));
    end
    if nargout ~= 1
        error('%s: Must have 1 output when using 2 inputs.', mfilename);
    end
    varargout{1} = [varargin{1}(:) num2cell(varargin{2})];
else
    error('%s: Check the inputs.', mfilename);
end