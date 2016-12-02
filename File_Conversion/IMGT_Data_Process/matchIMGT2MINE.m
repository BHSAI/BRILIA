%matchIMGT2MINE will take IMGT's VDJdata and MINE VDJdata, and then return
%only the resolved ones. This is so we can do fair Apple to Apple
%comparison on the QUALITY of annotation.


function matchIMGT2MINE()
%Determine the file sources
[FileName1, FilePath1] = uigetfile('*.xlsx','Select the IMGT format','*.IMGT');
[FileName2, FilePath2] = uigetfile('*.xlsx','Select the Mine format','*.T123');

[VDJdata1, NewHeader, ~, ~] = openSeqData([FilePath1 FileName1]);
[VDJdata2, ~, ~, ~] = openSeqData([FilePath2 FileName2]);

%Determine the file output names
DotLoc1 = find(FileName1 == '.');
DotLoc1 = DotLoc1(end);
SaveName1 = [FileName1(1:DotLoc1) ''];

DotLoc2 = find(FileName2 == '.');
DotLoc2 = DotLoc2(end);
SaveName2 = [FileName2(1:DotLoc2) ''];

%Select only the fully resolved VDJ segments of IMGT.
SeqNumLoc = findHeader(NewHeader,'SeqNum');
GrpNumLoc = findHeader(NewHeader,'GroupNum');
LengthLoc = findHeader(NewHeader,{'vAlignLength' 'dAlignLength' 'jAlignLength'});
FamNumLoc = findHeader(NewHeader,{'vMapNum','dMapNum','jMapNum'});



VDJresolved1 = cell2mat(VDJdata1(:,LengthLoc));
VDJresolved2 = cell2mat(VDJdata2(:,LengthLoc));

ValidLoc1 = ~isnan(sum(VDJresolved1,2)) & ~isempty(sum(VDJresolved1,2));
ValidLoc2 = ~isnan(sum(VDJresolved2,2)) & ~isempty(sum(VDJresolved2,2));
VDJdata1 = VDJdata1(ValidLoc1,:);
VDJdata2 = VDJdata2(ValidLoc2,:);

SeqNums1 = cell2mat(VDJdata1(:,SeqNumLoc));
SeqNums2 = cell2mat(VDJdata2(:,SeqNumLoc));
[~, IDX1, IDX2] = intersect(SeqNums1,SeqNums2);

VDJdata1 = VDJdata1(IDX1,:);
VDJdata2 = VDJdata2(IDX2,:);

for d1 = 1:size(VDJdata1,1)
    for d2 = 1:3
        VDJdata1{d1,FamNumLoc(d2)} = mat2str(VDJdata1{d1,FamNumLoc(d2)});
    end
end
for d1 = 1:size(VDJdata2,1)
    for d2 = 1:3
        VDJdata2{d1,FamNumLoc(d2)} = mat2str(VDJdata2{d1,FamNumLoc(d2)});
    end
end

%Save the extracted data as a file.
if ispc
    xlswrite([FilePath1 SaveName1 'IMGTnMINE.xlsx'],cat(1,NewHeader,VDJdata1));
    xlswrite([FilePath2 SaveName2 'IMGTnMINE.xlsx'],cat(1,NewHeader,VDJdata2));
else
    writeDlmFile(cat(1,NewHeader,VDJdata1),[FilePath1 SaveName1 'IMGTnMINE.csv'],'\t');
    writeDlmFile(cat(1,NewHeader,VDJdata2),[FilePath2 SaveName2 'IMGTnMINE.csv'],'\t');
end