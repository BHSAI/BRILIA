[FileNames, FilePath] = uigetfile('*Fix.xlsx','Select the file','multiselect','on');

for j = 1:length(FileNames)
    [~,~,FileData] = xlsread([FilePath, FileNames{j}]);
    if j == 1
        AllData = cell(100000,size(FileData,2));
        AllData(1:size(FileData,1),:) = FileData;
        q = size(FileData,1)+1;
    else
        AllData(q:q+size(FileData,1)-2,:) = FileData(2:end,:);
        q = q+size(FileData,1)-1; %removing headers
    end
end
AllData(q:end,:) = [];

[SaveName, SavePath] = uiputfile('*.xlsx','Save as');
xlswrite([SavePath SaveName],AllData);
