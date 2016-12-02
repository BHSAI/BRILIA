%Renaming Files
FileList = dir('*.txt');
for j = 1:size(FileList)
    FileName = FileList(j).name;
    NewFileName = strrep(FileName,'txt','csv');
    if strcmpi(FileName,NewFileName) == 0
        movefile(FileName,NewFileName);    
    end
end
