%removeDupSeq will look for duplicate sequences, combine them, and remove
%the extra entries. Duplicate sequences can arise after Indel correction or
%after padding/trimming sequences of variable lengths.
%
%  VDJdata = removeDupSeq(VDJdata, VDJheader)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    VDJheader: main BRILIA header cell

function VDJdata = removeDupSeq(VDJdata, VDJheader)
[H, L, Chain] = getAllHeaderVar(VDJheader);
GrpNum = cell2mat(VDJdata(:, H.GrpNumLoc));
UnqGrpNum = unique(GrpNum);
DelIdx = zeros(size(VDJdata, 1), 1, 'logical');
for y = 1:length(UnqGrpNum)
    if ~mod(y, 100)
        showStatus(sprintf('  Removing duplicate seq %d / %d.', y, length(UnqGrpNum)), []);
    end
    GrpIdx = GrpNum == UnqGrpNum(y);
    Tdata = VDJdata(GrpIdx, :);
    AncMapS = getTreeData(Tdata, VDJheader);
    AncMap = AncMapS.HAM;
    if any(AncMap(:, 3) == 0)
        Tdata = removeDupSeqPerGroup(VDJdata(GrpIdx, :), H, L, Chain);
        if sum(GrpIdx) ~= size(Tdata, 1)
            GrpLoc = find(GrpIdx);
            DelNum = sum(GrpIdx) - size(Tdata, 1);
            DelIdx(GrpLoc(end-DelNum+1:end)) = 1;
            VDJdata(GrpLoc(1:end-DelNum), :) = Tdata;
        end
    end
end
if max(DelIdx) > 0
    VDJdata(DelIdx, :) = [];
    showStatus(  'Deleted %d duplicate sequences.', sum(DelIdx));
end