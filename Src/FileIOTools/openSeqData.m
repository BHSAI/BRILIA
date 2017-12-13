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

if nargin > 0 && ischar(varargin{1}) && exist(varargin{1}, 'file')
    FullFileName = varargin{1};
elseif nargin == 0 || isempty(varargin{1})
    FullFileName = getBriliaFiles('', 'single');
    FullFileName = FullFileName{1};
else
    FullFileName = [];
end

if isempty(FullFileName)
    warning('%s: No valid file was choosen.', mfilename);
    return
end

[FilePath, FileName, FileExt] = parseFileName(FullFileName);
if isempty(FileExt) || isempty(FilePath)
    warning('%s: No valid file was found at "%s".', mfilename, FullFileName);
    return
end

[VDJheader, Delimiter] = readDlmFile(FullFileName, 'LineRange', [1 1]);
VDJinfo = readDlmFile('DataHeaderInfo.csv');
[~, Idx1, Idx2] = intersect(VDJheader, VDJinfo(:,1), 'stable');

if length(Idx1) == length(VDJheader) %You have all the info!
    IntLoc = find(startsWith(VDJinfo(Idx2, 2), 'int', 'ignorecase', true));
    StrPat = repmat('%s', 1, length(VDJheader));
    StrPat(IntLoc*2) = 'f';
    
    FID = fopen(FullFileName, 'r');
    VDJdata = textscan(FID, StrPat, 'Headerlines', 1, 'Delimiter', Delimiter, 'EndOfLine', '\r\n', 'CollectOutput', false);
    fclose(FID);
    VDJdata(:, IntLoc) = cellfun(@num2cell, VDJdata(:, IntLoc), 'unif', false);
    VDJdata = [VDJdata{:}];
    
    %Attempt to convert numeric matrix strings to numbers
    MatLoc = startsWith(VDJinfo(Idx2, 2), 'mat', 'ignorecase', true);
    VDJdata(:, MatLoc) = cellfun(@convStr2Num, VDJdata(:, MatLoc), 'unif', false);
else
    StrPat = repmat('%s', 1, length(VDJheader));
    
    FID = fopen(FullFileName, 'r');
    VDJdata = textscan(FID, StrPat, 'Headerlines', 1, 'Delimiter', Delimiter, 'EndOfLine', '\r\n', 'CollectOutput', true);
    fclose(FID);
    VDJdata = VDJdata{1};
    
    %Attempt to convert all numeric strings to numbers
    [~, ~, ~, NumLoc, ~] = getAllHeaderVar(VDJheader);
    for j = 1:length(NumLoc)
        for r = 1:size(VDJdata, 1)
            if isnumeric(VDJdata{r, NumLoc(j)}); continue; end
            try 
                VDJdata{r, NumLoc(j)} = convStr2Num(VDJdata{r, NumLoc(j)});
            catch 
            end
        end
    end
end

if nargout >= 5
    Map = getVDJmapper(VDJheader);
end