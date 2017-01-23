%readDlmFile will take a delimited file and parse it into a cell array of
%strings.
%
%  CellData = readDlmFile()
%
%  CellData = readDlmFile('Delimiter',Delimiter)
%
%  CellData = readDlmFile(InputFile)
%
%  CellData = readDlmFile(InputFile,'Delimiter',Delimiter)
%
%  INPUT
%    InputFile: full name of the file to be open
%    Delimiter [';' ',' '\t]: Separater marker between data. Default ';'.
%
%  OUTPUT
%    CellData: cell array of data
%
%  See also writeDlmFile

function CellData = readDlmFile(varargin)
%Look for delimeter input
Delimiter = ';';
for j = 1:length(varargin)
    if strcmpi(varargin{j},'delimiter');
        Delimiter = varargin{j+1};
        break
    end
end
varargin(j:j+1) = [];

%Get the input file name
if isempty(varargin)
    [FileName, FilePath] = uigetfile('*.csv','Open CSV file');
    InputFileName = [FilePath FileName];
else
    if ~strcmpi(varargin{1},'delimiter');
        InputFileName = varargin{1};
    else
        [FileName, FilePath] = uigetfile('*.csv','Open CSV file');
        InputFileName = [FilePath FileName];
    end
end

%Determine line count first
FID = fopen(InputFileName,'r');
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
Header = regexp(HeaderTxt,Delimiter,'split');

%Start the matrix generation
CellData = cell(NumLines,length(Header));
FID = fopen(InputFileName,'r');
for j = 1:NumLines
    TextLine = fgetl(FID);
    TextLine = strrep(TextLine,'"',''); %For some reason, fgetl will return a ' " ' at the beginning and end of line.
    if ~strcmpi(Delimiter,',')
        TextLine = strrep(TextLine,',',''); %For some reason, there can be string of ",,,,,," at end of file.
    end
    LineData = regexp(TextLine,Delimiter,'split');
    CellData(j,1:length(LineData)) = LineData;
end
fclose(FID);
