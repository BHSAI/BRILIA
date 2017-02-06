%extractGeneFamily will open up adaptive file formats, and extract the
%family number, placing them into a csv file with " ; " delimiter, of V; D;
%J family # only. EX: IGHV01, IGHD02r, IGHJ01 will be 1; -2; 1.

function extractGeneFamily()

[FileName, FilePath] = uigetfile('*.xlsx','Open Excel file Adap format');

[~, ~, SampleData] = xlsread([FilePath FileName]);
[SampleData, Header] = filterHeader(SampleData);

%Look for the V,D,J family name column in the spreadsheet
VDJfamNums = zeros(size(SampleData,1),3);
VDJcolNums = zeros(1,3);
for j = 1:length(Header)
    if strcmpi(Header{j},'vFamilyName')
        VDJcolNums(1) = j;
    end
    
    if strcmpi(Header{j},'dFamilyName')
        VDJcolNums(2) = j;
    end
    
    if strcmpi(Header{j},'jFamilyName')
        VDJcolNums(3) = j;
    end
    
end

for k = 1:length(VDJcolNums);
    for j = 1:size(SampleData,1)
        NameQuery = SampleData{j,VDJcolNums(k)};
        if isnan(NameQuery) 
            FamNum = 0;
        elseif isempty(NameQuery)
            FamNum = 0;
        else
            FamNum = str2double(regexp(SampleData{j,VDJcolNums(k)},'\d*','match'));
            if isempty(FamNum) %In case it's "Unresolved"
                FamNum = 0;
            end
            if isempty(regexp(SampleData{j,VDJcolNums(k)},'r','once')) == 0
                FamNum = FamNum*-1;
            end
        end
        
        VDJfamNums(j,k) = FamNum;
    end
end

%Save the file as a csv
DotLoc = regexp(FileName,'\.');
DotLoc = DotLoc(end);
OutputFileName = [FileName(1:DotLoc-1) '.csv'];
dlmwrite([FilePath OutputFileName],VDJfamNums,'delimiter',';');

