%readDlmFile will take a delimited file of anytime and parse it into a cell
%array of strings.
%
%
function CellData = readDlmFile(varargin)
if isempty(varargin)
    [FileName, FilePath] = uigetfile('*.csv','Open CSV file');
    FullName = [FilePath FileName];
else
    if ~strcmpi(varargin{1},'delimiter');
        FullName = varargin{1};
    else
        [FileName, FilePath] = uigetfile('*.csv','Open CSV file');
        FullName = [FilePath FileName];
    end
end

%Look for delimeters
D = '\t'; %by default.
for j = 1:length(varargin)
    if strcmpi(varargin{j},'delimiter');
        D = varargin{j+1};
    end
end

%Determine line count first
FID = fopen(FullName,'r');
NumLines = 0;
HeaderTxt = '';
while 1
    TextLine= fgetl(FID);
    if ischar(TextLine)
        if isempty(HeaderTxt)
            HeaderTxt = TextLine;
        end
        NumLines = NumLines+1;
    else
        break
    end
end
fclose(FID);

%Parse the headers (just to get the column count);
Header = regexp(HeaderTxt,D,'split');

%Start the matrix generation
CellData = cell(NumLines,length(Header));
FID = fopen(FullName,'r');
for j = 1:NumLines
    TextLine = fgetl(FID);
    TextLine = strrep(TextLine,'"',''); %For some reason, fgetl will return a ' " ' at the beginning and end of line.
    if ~strcmpi(D,',')
        TextLine = strrep(TextLine,',',''); %For some reason, there can be string of ",,,,,," at end of file.
    end
    LineData = regexp(TextLine,D,'split');
    CellData(j,1:length(LineData)) = LineData;
end
fclose(FID);
