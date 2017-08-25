%readDlmFile will take a delimited file and parse it into a cell array of
%strings.
%
%  CellData = readDlmFile()
%
%  CellData = readDlmFile(InputFile)
%
%  CellData = readDlmFile(..., 'Delimiter', Delimiter)
%
%  CellData = readDlmFile(..., 'LineRange', LineRange)
%
%  INPUT
%    InputFile: full name of the file to be open
%    Delimiter ['' ';' ',' '\t']: Delimiter for the file. The default is
%      empty '', meaning it will autodetect based on first two text lines.
%      The most common delimiter symbol is assumed the true delimiter.
%    LineRange [1 Inf]: 1x2 matrix of start to end lines to extract
%
%  OUTPUT
%    CellData: cell array of string data
%
%  WARNING
%    If the delimited file has improper delimiter usage, such as a single
%    string field having a delimiter symbol, then will split it
%    incorrectly. Also, if the number of delimited column is different
%    between data rows, this will delete off the excess data entries,
%    mainly because these are often OS-specific addons of ',' and to
%    prevent indexing errors and rescaling of the intended data table.
%
%  See also writeDlmFile

function CellData = readDlmFile(varargin)
CellData = {};

%Look for Delimiter and LineRange inputs
Delimiter = '';
LineRange = [1 Inf];
SearchFor = ''; %regexp pattern to look for
SearchAt = []; %the nth column to look in
DelThis = zeros(1, length(varargin), 'logical');
for j = 1:length(varargin)
    if ischar(varargin{j})
        if strcmpi(varargin{j}, 'Delimiter')
            Delimiter = varargin{j+1};
            if ~ischar(Delimiter)
                error('%s: Delimiter not valid', mfilename);
            end
            DelThis(j:j+1) = 1;
        end
        if strcmpi(varargin{j}, 'LineRange')
            LineRange = varargin{j+1};
            if ~isnumeric(LineRange)
                error('%s: LineRange not valid', mfilename);
            end
            DelThis(j:j+1) = 1;
        end
        if strcmpi(varargin{j}, 'SearchFor')
            SearchFor = varargin{j+1};
            if ~ischar(SearchFor)
                error('%s: SearchFor must be a string.', mfilename);
            end
            DelThis(j:j+1) = 1;
        end
        if strcmpi(varargin{j}, 'SearchAt')
            SearchAt = varargin{j+1};
            if ~isnumeric(SearchAt)
                error('%s: SearchAt not valid', mfilename);
            end
            DelThis(j:j+1) = 1;
        end
    end
end
varargin(DelThis) = [];

%Get the input file name
if isempty(varargin) || (~isempty(varargin) && (~ischar(varargin{1})))
    [FileName, FilePath] = uigetfile('*.csv;*.tsv;*.ssv', 'Open delimited file.');
    InputFileName = [FilePath FileName];
else
    InputFileName = varargin{1};
end
if ~exist(InputFileName, 'file')
    return;
end

%Determine line count and the first two txt
[FID, MSG] = fopen(InputFileName, 'r');
if FID < 0
    error('%s: Could not open file [ %s ].\n  %s', mfilename, InputFileName, MSG);
end
NumLines = 0;
TextLine = fgetl(FID);
Line1Txt = '';
if ischar(TextLine)
    Line1Txt = TextLine;
    NumLines = NumLines + 1;
end
TextLine = fgetl(FID);
Line2Txt = '';
if ischar(TextLine)
    Line2Txt = TextLine;
    NumLines = NumLines + 1;
end
while ~feof(FID)
    fgetl(FID);
    NumLines = NumLines+1;
end

%Determine the delimiter if unknown
if isempty(Delimiter)
    DelimiterChoice = {',' ';' '\t'};
    Line1Ct = zeros(1, length(DelimiterChoice));
    Line1Ct(1) = length(regexp(Line1Txt, '\,'));
    Line1Ct(2) = length(regexp(Line1Txt, '\;'));
    Line1Ct(3) = length(regexp(Line1Txt, '\s+'));
    if isempty(Line2Txt)
        Line2Ct = Line1Ct;
    else
        Line2Ct = zeros(1, length(DelimiterChoice));
        Line2Ct(1) = length(regexp(Line2Txt, '\,'));
        Line2Ct(2) = length(regexp(Line2Txt, '\;'));
        Line2Ct(3) = length(regexp(Line2Txt, '\s+'));
    end
    
    %Delimiter is the most frequent one with equal number in both lines
    SortMat = [Line1Ct' Line2Ct' [1:3]'];
    ValidOnes = Line1Ct' == Line2Ct';
    SortMat = SortMat(ValidOnes, :);
    SortMat = sortrows(SortMat);
    Delimiter = DelimiterChoice{SortMat(end, 3)};
end

%Seek the correct line specified by the LineRange
fseek(FID, 0, 'bof');
if LineRange(1) < 1; LineRange = 1; end
if LineRange(2) > NumLines; LineRange(2) = NumLines; end
for k = 1:LineRange(1)-1
    fgetl(FID);
end

%Extract the requested data
ColCt = length(regexp(Line1Txt, Delimiter, 'split'));
CellData = cell(diff(LineRange)+1, ColCt);
if isempty(SearchFor)
    for j = 1:diff(LineRange)+1
        TextLine = fgetl(FID);
        LineData = regexp(TextLine, Delimiter, 'split');
        MinColCt = min(length(LineData), ColCt);
        CellData(j, 1:MinColCt) = LineData(1:MinColCt);
    end
else
    if isempty(SearchAt) %Search everywhere
        SearchAt = 1:ColCt;
    end
    j = 1;
    for k = 1:diff(LineRange)+1
        TextLine = fgetl(FID);
        LineData = regexp(TextLine, Delimiter, 'split');
        KeepThis = 0;
        for q = 1:length(SearchAt)
            if ~isempty(regexpi(LineData{SearchAt(q)}, SearchFor, 'once'))
                KeepThis = 1;
                break;
            end
        end
        if KeepThis
            MinColCt = min(length(LineData), ColCt);
            CellData(j, 1:MinColCt) = LineData(1:MinColCt);
            j = j + 1;
        end
    end
    CellData(j:end, :) = [];
end
fclose(FID);
