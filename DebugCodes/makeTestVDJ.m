%Examples
[VDJdata,NewHeader] = convertInput2VDJdata([cd '\TestFullVDJ_RepSize1000.xlsx']);
getHeaderVar;
for j = 1:30 %Constant region of J
    VDJdata{j,SeqLoc} = [VDJdata{j,SeqLoc} 'ACGCTG'];
end
for j = 31:60 %Barcode of V
    VDJdata{j,SeqLoc} = ['ACGCGT' VDJdata{j,SeqLoc}];
end
for j = 61:90 %Random character in seq
    VDJdata{j,SeqLoc}(3) = '-';
end
for j = 91:120 %Rev complement
    VDJdata{j,SeqLoc} = seqrcomplement(VDJdata{j,SeqLoc});
end
for j = 121:150 %Add nts here and there
    VDJdata{j,SeqLoc} = seqrcomplement(VDJdata{j,SeqLoc});
end
VDJdata(41:end,:) = [];

writeDlmFile([NewHeader;VDJdata],[[cd '\'] 'FullVDJtest.csv'],'\t');
clear
clc