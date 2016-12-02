%standardizeAdapFiles will take all files associated with Adap biotech from
%Chris Coopers data, and then convert the files to be of standard format.
%This is needed to simplify data analysis.

function standardizeAdapFiles()
disp('WARNING: This will process all tsv files');
while 1
    ContThis = input('Continue? y or n','s');
    if strcmpi(ContThis,'y');
        break
    elseif strcmpi(ContThis,'n');
        break
    end
end
if strcmpi(ContThis,'n'); return; end

%Renaming tsv files as csv files for MatLab
FileList = dir('*.tsv*');
for j = 1:size(FileList)
    FileName = FileList(j).name;
    NewFileName = strrep(FileName,'.tsv','.csv');
    if strcmpi(FileName,NewFileName) == 0
        movefile(FileName,NewFileName);    
    end
end

%Rearranging columns tsv files as csv files for MatLab
[~, ~, StandardData] = xlsread('Headers.xlsx');
AdapHeaderLoc = findHeader(StandardData(1,:),'Adap');
StandardHeader = StandardData(2:end,AdapHeaderLoc);

FileList = dir('*.csv*');
for j = 1:size(FileList)
    FileName = FileList(j).name;
    SampleData = readDlmFile(FileName,'delimiter','\t');
    [SampleData, SampleHeader] = filterHeader(SampleData);
    
    Map1 = 1:length(SampleHeader);
    Map2 = findHeader(StandardHeader,SampleHeader);
    DelMap = Map2 == 0;
    Map1(DelMap) = [];
    Map2(DelMap) = [];
    
    %Generate the cell matrix and reshuffle the table.
    NewTable = cell(size(SampleData,1),length(StandardHeader));
    NewTable(:,Map2) = SampleData(:,Map1);
    NewTable = cat(1,StandardHeader',NewTable);
       
    %Outputfile
    DotLoc = regexp(FileName,'\.');
    FilePre = FileName(1:DotLoc(end)-1);
    writeDlmFile(NewTable,[FilePre '.std.csv'],'\t');
end

%Move the original files into Original folder
if exist('Original','dir') == 0
    mkdir('Original');
end
for k = 1:length(FileList)
    movefile(FileList(k).name,'Original')
end