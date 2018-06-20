%combineDlmFile will simply combine delimited files together, but only if
%there are the same number of columns in each file.
%
%  combineDlmFile(FileList, OutputName, )
%
%  INPUT
%    FileList: Mx1 cell of file names to combine
%    OutputName: the file name to store the combined files
%  
%  OUTPUT
%    Will generated a combined delimited file, separated by a single line
%    spacer. A -1 would mean to exclude the 1st line of all subsequent
%    files being combined, which is good if you want to avoid the header
%    lines.
function combineDlmFile(FileList, OutputName, SpacerNum)
if nargin < 1 || isempty(FileList)
    FileList = openFileDialog('*.*sv', 'multiselect', 'on');
end

if nargin < 2 || isempty(OutputName)
    [OutFileName, OutFilePath] = uiputfile('CombinedFile.csv', 'Save combined file as');
    OutputName = fullfile(OutFilePath, OutFileName);
end

if nargin < 3 || isempty(SpacerNum)
    SpacerNum = 0;
end

%Make sure there's equal number of columns
ColNum = 0;
for j = 1:length(FileList)
    [ColNames, Delimiter] = readDlmFile(FileList{j}, 'LineRange', [1 1]);
    if j == 1
        ColNum = length(ColNames);
    elseif ColNum ~= length(ColNames)
        error('%s: All file must have the same number of columns to append to one file.', mfilename);
    end
end

if SpacerNum >= 0
    SpacerLine = repmat({''}, SpacerNum, ColNum);
    LineRange = [1 Inf];
else
    SpacerLine = {};
    LineRange = [abs(SpacerNum) Inf];
end

%Proceed to combine
for j = 1:length(FileList)
    Data = readDlmFile(FileList{j}, 'LineRange', LineRange);
    if j == 1
        writeDlmFile(Data, OutputName, Delimiter);
    else
        writeDlmFile([SpacerLine; Data], OutputName, Delimiter, 'append');
    end
end