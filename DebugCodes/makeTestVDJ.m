%Examples
[VDJdata,VDJheader] = convertInput2VDJdata([cd '\TestFullVDJ_RepSize1000.xlsx']);
H = getHeaderVar(VDJheader);
for j = 1:30 %Constant region of J
    VDJdata{j,H.SeqLoc} = [VDJdata{j,H.SeqLoc} 'ACGCTG'];
end
for j = 31:60 %Barcode of V
    VDJdata{j,H.SeqLoc} = ['ACGCGT' VDJdata{j,H.SeqLoc}];
end
for j = 61:90 %Random character in seq
    VDJdata{j,H.SeqLoc}(3) = '-';
end
for j = 91:120 %Rev complement
    VDJdata{j,H.SeqLoc} = seqrcomplement(VDJdata{j,H.SeqLoc});
end
for j = 121:150 %Add nts here and there
    VDJdata{j,H.SeqLoc} = seqrcomplement(VDJdata{j,H.SeqLoc});
end
VDJdata(41:end,:) = [];

writeDlmFile([VDJheader;VDJdata],[[cd '\'] 'FullVDJtest.csv'],'\t');
clear
clc
