%openSeqData will ask user to open up a csv or xlsx file, and then return
%the data and header information.

function [SampleData, Header, FileName, FilePath] = openSeqData(varargin)
%Look for the file and parse file name
if isempty(varargin)
    [FileName, FilePath] = uigetfile('*.xlsx;*.csv;*.fa;*.fastq','Find SampleData format file');
    FullName = [FilePath FileName];
else
    FullName = varargin{1};
    if isempty(FullName)
        [FileName, FilePath] = uigetfile('*.xlsx;*.csv;*.fa;*.fastq','Find SampleData format file');
        FullName = [FilePath FileName];
    end
end
[FilePath, FileName, FileExt] = parseFileName(FullName);

%Determine if  you want to eval everything.
EvalThis = 'eval';
if length(varargin) == 2
    EvalThis = varargin{2};
end

%Load the file
if strcmpi(FileExt,'.csv');
    SampleData = readDlmFile(FullName,'delimiter',';');
    if size(SampleData,2) == 1 %Probably wrong delimiter
        SampleData = readDlmFile(FullName,'delimiter','\t');
    end
elseif strcmpi(FileExt,'.xlsx')
    [~, ~, SampleData] = xlsread(FullName);
else
    error('Unrecognized file extension.');
end

%Split the file into the Data table and Header table.
if length(varargin) == 3
    HeaderCt = varargin{3};
else
    HeaderCt = 1;
end
[SampleData, Header] = filterHeader(SampleData,HeaderCt);

%Attempt to convert all evaluable strings to numbers
if strcmpi(EvalThis,'eval') %Evaluate all cells that are matlab variables.
    for r = 1:size(SampleData,1)
        for c = 1:size(SampleData,2)
            if ischar(SampleData{r,c})
                if isempty(SampleData{r,c}); continue; end
                if isnan(SampleData{r,c}); continue; end
                if isempty(regexp(SampleData{r,c},'[^\d*\]\[\s]','once'))
                    SampleData{r,c} = eval(SampleData{r,c});
                end
            end
        end
    end
end
