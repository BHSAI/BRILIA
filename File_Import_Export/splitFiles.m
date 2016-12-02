[FileName, FilePath] = uigetfile('*.xlsx','Select the file');
[~,~,FileData] = xlsread([FilePath, FileName]);
[FileData, Header] = filterHeader(FileData,1);

SeqFileLoc = findHeader(Header,'SeqFile');

[UnqFile, ~, UnqFileIdx] = unique(FileData(:,SeqFileLoc));
for j = 1:length(UnqFile)
    xlswrite([FilePath UnqFile{j} '.xlsx'],[Header;FileData(UnqFileIdx == j,:)]);
end
