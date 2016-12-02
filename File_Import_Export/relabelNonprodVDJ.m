%relabelNonprodVDJ simply performs and saves the VDJdata files using the
%current labeling procedure for productive, nonproductive, or maybe
%productive sequences, which is Y, N, and M respectively in the funcational
%column.

%open and label VDJdata.
[VDJdata,NewHeader,FileName,FilePath] = openSeqData();
getHeaderVar;
VDJdata = labelNonprodVDJ(VDJdata,NewHeader);

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
    xlswrite([FilePath SaveName '.xlsx'],[NewHeader; VDJdata]);
else
    writeDlmFile([NewHeader;VDJdata],[FilePath SaveName '.csv'],'\t');
end
