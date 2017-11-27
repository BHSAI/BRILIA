%splitFileByGroup will open a dlm file and split it based on the unique
%column or conatenated columns. 
%
%  INPUT
%    FileName: filename of the dlm file
%    ColNum: column numbers to use to determine unique groups
%
%  OUTPUT
%    Output are files with the same filename, but with a group number
%    assigned to it starting from 1 to N. Example filename_1.csv,
%    filename_2.csv, filename_3.csv,  etc...

function splitFileByGroup(FileName, ColNum)
if ~exist(FileName, 'file') && ~exist(fullfile(pwd, FileName), 'file')
    error('%s: Could not find file "%s".', mfilename, FileName);
end

Data = readDlmFile(FileName);
Header = Data(1, :);
Data(1, :)  = [];
[UnqName, ~, UnqIdx] = unique(Data(:, ColNum));

[FilePath, InFileName, FileExt] = parseFileName(FileName);
DotLoc = find(InFileName == '.');
OutPre = InFileName(1:DotLoc(end)-1);

for k = 1:max(UnqIdx)
    Idx = UnqIdx == k;
    OutFileName = fullfile(FilePath, [OutPre '_' UnqName{k} FileExt]);
    writeDlmFile(cat(1, Header, Data(Idx, :)), OutFileName, ',');
end