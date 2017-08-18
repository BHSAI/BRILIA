%conformGeneGroup will force all sequences in the same groups to have the
%same Gene family and VMDNJ segment lengths as the first sequence in the
%cluster, which is suppose to be close to the germline sequence.
%
%  VDJdata = conformGeneGroup(VDJdata, VDJheader, DB)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    VDJheader: main BRILIA header cell
%    DB: Gene database structure (getCurrentDatabase.m)
%
%  OUTPUT
%    VDJdata: modified VDJdata where all seq in a cluster has the same
%      gene annotations

function VDJdata = conformGeneGroup(VDJdata, VDJheader, DB)
[H, L, Chain] = getAllHeaderVar(VDJheader);

%Standadize annotations
GrpNum = cell2mat(VDJdata(:, H.GrpNumLoc));
UnqGrpNum = unique(GrpNum);
for y = 1:length(UnqGrpNum)
    IdxLoc = UnqGrpNum(y) == GrpNum;
    if sum(IdxLoc) <= 1; continue; end %Nothing to correct, single entry.
    Tdata = VDJdata(IdxLoc, :);    
    
    if H.SeqLoc > 0 %Heavy chain
        %Standardize the other fields as well (1st one is the root)
        Tdata(:, [H.LengthLoc(:); H.GeneNameLoc(:); H.GeneNumLoc(:); H.DelLoc(:)]) = repmat(Tdata(1, [H.LengthLoc(:); H.GeneNameLoc(:); H.GeneNumLoc(:); H.DelLoc(:)]), size(Tdata, 1), 1);
    end
    
    if L.SeqLoc > 0 %Light chain
        %Standardize the other fields as well (1st one is the root)
        Tdata(:, [L.LengthLoc(:); L.GeneNameLoc(:); L.GeneNumLoc(:); L.DelLoc(:)]) = repmat(Tdata(1, [L.LengthLoc(:); L.GeneNameLoc(:); L.GeneNumLoc(:); L.DelLoc(:)]), size(Tdata, 1), 1);
    end
    
    VDJdata(IdxLoc, :) = Tdata;  
end

%Update all data
VDJdata = buildRefSeq(VDJdata, VDJheader, DB, Chain, 'germline', 'first'); %must do first seq of each cluster
VDJdata = updateVDJdata(VDJdata, VDJheader, DB);
