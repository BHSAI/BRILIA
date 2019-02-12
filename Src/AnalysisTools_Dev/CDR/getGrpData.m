%getGrpData will collect the size of the clonotypes from VDJdata.
%
%  Data = getGrpData(VDJdata, VDJheader)
%
%  INPUT
%    VDJdata: BRILIA data cell
%    VDJheader: BRILIA header cell
%
%  OUTPUT
%    Data: structure of group size data with fields:
%      GroupNum - group number
%      Template - sum of template counts of that group
%      Members  - number of unique B cell clones

function Data = getGrpData(VDJdata, VDJheader)
Map = getVDJmapper(VDJheader);
GrpNum = cell2mat(VDJdata(:, Map.GrpNum));
UnqGrpNum = unique(GrpNum);

Data(length(UnqGrpNum)) = struct('GroupNum', 0, 'Template', 0, 'Members', 0, 'HCDR3', '', 'LCDR3', '', 'HCDR3Count', 0, 'LCDR3Count', 0);
for y = 1:length(UnqGrpNum)
    Idx = find(GrpNum == UnqGrpNum(y));
    Data(y).GroupNum = VDJdata{Idx(1), Map.GrpNum};
    Data(y).Template = sum(cell2mat(VDJdata(Idx, Map.Template)));
    Data(y).Members = length(Idx);
    if contains(Map.Chain, 'H')
        Data(y).HCDR3 = unique(VDJdata(Idx, Map.hCDR3(1)));
        Data(y).HCDR3Count = length(Data(y).HCDR3);
    end
    if contains(Map.Chain, 'L')
        Data(y).UnqLCDR3 = unique(VDJdata(Idx, Map.lCDR3(1)));
        Data(y).LCDR3Count = length(Data(y).LCDR3);
    end
end