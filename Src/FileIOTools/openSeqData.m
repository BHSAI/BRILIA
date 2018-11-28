%openSeqData will ask user to open up a delimited BRILIA output file for
%use in matlab. Determines which data should be numerical and converts
%strings to numbers accordingly.
% 
%  [VDJdata, VDJheader] = openSeqData
%
%  [VDJdata, VDJheader] = openSeqData(FileName)
%
%  [VDJdata, VDJheader, FileName, FilePath] = openSeqData(...)
%
%  [VDJdata, VDJheader, FileName, FilePath, Map] = openSeqData(...)
%
%  INPUT
%    FileName: the .csv file from the BRILIA output
%
%  OUTPUT
%    VDJdata: main BRILIA annotation cell array
%    VDJheader: header anmes of VDJdata
%    FileName: the output file name (used if file is selected by GUI)
%    FilePath: the output file path (used if file is selected by GUI)
%    Map: a structure relating VDJheader to column number (used for coding)
function [VDJdata, VDJheader, FileName, FilePath, Map] = openSeqData(varargin)
VDJdata = [];
VDJheader = [];
FileName = '';
FilePath = '';
Map = [];

if nargin == 0
    FullFileName = getBriliaFiles('', 0);
else
    if contains(varargin{1}, {'.', filesep})
        FullFileName = dir2(varargin{1}, 'file');
        varargin = varargin(2:end);
    else
        FullFileName = getBriliaFiles('', 'single');
    end
end
if isempty(FullFileName)
    warning('%s: No valid file was choosen.', mfilename);
    return
elseif numel(FullFileName) > 1
    warning('%s: Cannot select multiple files. Choosing 1st one "%s".', mfilename, FullFileName{1});
end
FullFileName = FullFileName{1}; %Unwrap 1st one only

[FilePath, FileName, FileExt] = parseFileName(FullFileName);
if isempty(FileExt) || isempty(FilePath)
    warning('%s: No valid file was found at "%s".', mfilename, FullFileName);
    return
end

[VDJheader, Delimiter] = readDlmFile(FullFileName, 'LineRange', 1);

%Determine numeric, matrix, and string columns
[Map, NumIdx] = getVDJmapper(VDJheader);
MatIdx = find(endsWith(VDJheader, 'MapNum', 'ignorecase', true)); %Determine the 'MapNum' columns, as these require string to matrix conversion. Will be deprecated in future.
NumIdx = setdiff(NumIdx, MatIdx);

%Setup the pattern search
StrPat = repelem({'%s'}, 1, length(VDJheader));
StrPat(NumIdx) = {'%f'};
StrPat = [StrPat{:}];

%Read in the data and format the columns into a single cell array
FID = fopen(FullFileName, 'r');
VDJdata = textscan(FID, StrPat, 'Headerlines', 1, 'Delimiter', Delimiter, 'EndOfLine', '\r\n', 'CollectOutput', false);
fclose(FID);
VDJdata(:, NumIdx) = cellfun(@num2cell, VDJdata(:, NumIdx), 'unif', false);
VDJdata(:, MatIdx) = cellfun(@convStr2NumMEX, VDJdata(:, MatIdx), 'unif', false);
VDJdata = [VDJdata{:}];

%Filter for relevant data
if ~(any(cellfun('isempty', varargin)))
    HeaderStr = formatStrSame(VDJheader);
    QueryStr = formatStrSame(varargin);
    [~, GetIdx] = intersect(HeaderStr, QueryStr);
    if ~isempty(GetIdx)
        VDJdata = VDJdata(:, GetIdx);
        VDJheader = VDJheader(:, GetIdx);
    end
end

%Format string to be same for matching purposes. 
%EDIT_NOTE: Edit this code if the string matching criteria changes
function Str = formatStrSame(Str)
try
    Str = lower(regexprep(Str, '[^a-zA-Z0-9\|]', ''));
catch
    Str = Str;
end