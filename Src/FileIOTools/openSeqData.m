%openSeqData will ask user to open up a delimited BRILIA output file for
%use in matlab. Determines which data should be numerical and converts
%strings to numbers accordingly.
% 
%  [VDJdata, VDJheader] = openSeqData;
%
%  [VDJdata, VDJheader, FileName, FilePath] = openSeqData;
%
%  ... = openSeqData(FileName, Param, Value, ...)
%
%  ... = openSeqData('', Param, Value, ...)

function [VDJdata, VDJheader, FileName, FilePath] = openSeqData(varargin)
%Default output for failed open
VDJdata = [];
VDJheader = [];
FileName = '';
FilePath = '';

%Look for the file and parse file name
if isempty(varargin) || (~isempty(varargin) && isempty(varargin{1}))
    [FileName, FilePath] = uigetfile('*.csv;*.tsv;*.ssv', 'Open BRILIA output file.');
    if isnumeric(FileName) %No file selected
        return;
    end
    FullName = [FilePath FileName];
    if ~isempty(varargin)
        varargin(1) = [];
    end
else
    FullName = varargin{1};
    if exist(FullName, 'file') %Check if this input was a file name
        varargin(1) = [];
    else %Probably was just a param-field
        [FileName, FilePath] = uigetfile('*.csv;*.tsv;*.ssv', 'Open BRILIA output file.');
        if isnumeric(FileName) %No file selected
            return;
        end
        FullName = [FilePath FileName];
    end
end
[FilePath, FileName, FileExt] = parseFileName(FullName);

%Parse the remaining inputs
P = inputParser;
addParameter(P, 'Delimiter', '', @(x) ismember({lower(x)}, {';', ',', '\t', ''}));
addParameter(P, 'SearchColName', '', @(x) ischar(x) || iscell(x));
addParameter(P, 'SearchFor', '', @(x) ischar(x)); 
parse(P, varargin{:});
P = P.Results;

%Read valid files
if isempty(FileExt) || isempty(FilePath)
    warning('%s: No valid file was found at [ %s ].', mfilename, FullName);
    return;
end
VDJheader = readDlmFile(FullName, 'Delimiter', P.Delimiter, 'LineRange', [1 1]);
if ~isempty(P.SearchFor)
    SearchAt = findCell(VDJheader, P.SearchColName, 'MatchCase', 'Any');
    if SearchAt(1) == 0
        SearchAt = 1:length(VDJheader);
    end
else
    SearchAt = [];
end
VDJdata = readDlmFile(FullName, 'Delimiter', P.Delimiter, 'LineRange', [2 Inf], 'SearchFor', P.SearchFor, 'SearchAt', SearchAt);

%Attempt to convert numeric strings to numbers
[~, ~, ~, NumLoc, ~] = getAllHeaderVar(VDJheader);
for j = 1:length(NumLoc)
    for r = 1:size(VDJdata, 1)
        try 
            VDJdata{r, NumLoc(j)} = eval(VDJdata{r, NumLoc(j)});
        catch 
        end
    end
end
