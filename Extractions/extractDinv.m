%extractDinv will open a file, extract the D inverse matches only, and save
%it separately. Useful for rerunning match using matchVDJ, but with D
%inverse turned off.

function extractDinv()

[VDJdata, VDJheader, FileName, FilePath] = openSeqData();

%Locate necessary data columns
H.FamNumLoc = findHeader(VDJheader,{'vMapNum','dMapNum','jMapNum'});

KeepThis = zeros(size(VDJdata,1),1) > 0;
for j = 1:size(VDJdata,1)
    FamNum = VDJdata{j,H.FamNumLoc(2)};
    for k = 1:length(FamNum)
        if mod(FamNum(k),2) == 0
            KeepThis(j) = 1;
            break;
        end
    end
end

VDJdata = VDJdata(KeepThis,:);

%Prepare to save the files into new folder
DotLoc = find(FileName == '.');
SaveNamePre = FileName(1:DotLoc(end)-1);
SavePath = FilePath;

%Before saving to xlsx, convert columns with matrix values into char for saving
VDJdataRaw = reformatAlignment(VDJdata,1);
for d1 = 1:size(VDJdata,1)
    for d2 = 1:3
        VDJdataRaw{d1,H.FamNumLoc(d2)} = mat2str(VDJdataRaw{d1,H.FamNumLoc(d2)});
    end
end

if ispc
    xlswrite([SavePath SaveNamePre '.Dinv.xlsx'],cat(1,VDJheader,VDJdataRaw));
else
    writeDlmFile(cat(1,VDJheader,VDJdataRaw),[SavePath SaveNamePre '.Dinv.csv'],'\t');
end
