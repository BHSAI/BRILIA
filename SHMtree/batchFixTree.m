%Open multiple files
[FileNames,FilePath] = uigetfile('*.xlsx;*.csv','Open files','multiselect','on');
if ischar(FileNames)
    FileNames = {FileNames};
end
for f = 1:length(FileNames)
    [VDJdata,NewHeader,FileName,FilePath] = openSeqData([FilePath FileNames{f}]);
    getHeaderVar;
    
    VDJdata = fixTree(VDJdata,NewHeader);
    
    %Save the files
    DotLoc = find(FileName == '.');
    DotLoc = DotLoc(end);
    SaveName = FileName(1:DotLoc-1);
    %Before saving to xlsx, convert columns with matrix values into char
    for q = 1:size(VDJdata,1)
        for w = 1:3
            VDJdata{q,FamNumLoc(w)} = mat2str(VDJdata{q,FamNumLoc(w)});
        end
    end
    if ispc
        xlswrite([FilePath SaveName 'D.xlsx'],[NewHeader; VDJdata]);
    else
        writeDlmFile([NewHeader;VDJdata],[FilePath SaveName 'D.csv'],'\t');
    end
    
end