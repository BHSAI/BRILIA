%This script will combine all the .Fix.xlsx files into one

%Renaming tsv files as csv files for MatLab
FileList = dir('*Fix.xlsx');
FilePath = [cd '\'];
FileName = FileList(1).name;

[VDJdata,NewHeader,~,~] = openSeqData([FilePath FileName]);
getHeaderVar;
FullData = cell(110000,length(NewHeader));
s1 = 1;
s2 = 0;

GrpNumCt = 1;
for j = 1:size(FileList,1)
    FileName = FileList(j).name;
    [VDJdata,NewHeader,~,~] = openSeqData([FilePath FileName]);
    s1 = s2+1;
    s2 = s1 + size(VDJdata,1) - 1;    
    
    %Renumber group number
    GrpNums = cell2mat(VDJdata(:,GrpNumLoc));
    GrpNum = zeros(size(GrpNums));
    [UnqGrp, ~, ~] = unique(GrpNums);
    for b = 1:length(UnqGrp)
        GrpLoc =  GrpNums == UnqGrp(b);
        GrpNum(GrpLoc) = GrpNumCt;
        GrpNumCt = GrpNumCt+1;
    end
    VDJdata(:,GrpNumLoc) = num2cell(GrpNum);    
    
    FullData(s1:s2,:) = VDJdata;
    j
end

FullData(s2+1:end,:) = [];

%Before saving to xlsx, convert columns with matrix values into char
FullData = reformatAlignment(FullData,1);
for q = 1:size(FullData,1)
    for w = 1:3
        FullData{q,FamNumLoc(w)} = mat2str(FullData{q,FamNumLoc(w)});
    end
    if mod(q,1000) == 0
        [q]
    end
end
% writeDlmFile([NewHeader;FullData],'FullList2.csv');
xlswrite('FullList2.xlsx',[NewHeader;FullData]);

%filterNonprodVDJ()