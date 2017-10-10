%spliceData will split VDJdata into smaller cells of VDJdata to enable
%parallel processing of groups.
function ParVDJdata = spliceData(VDJdata, VDJheader)
H = getAllHeaderVar(VDJheader);
GrpNum = cell2mat(VDJdata(:, H.GrpNumLoc));
UnqGrpNum = unique(GrpNum);
ParVDJdata = cell(1, length(UnqGrpNum));
for k = 1:length(UnqGrpNum)
    ParVDJdata{k} = VDJdata(UnqGrpNum(k) == GrpNum, :);
end