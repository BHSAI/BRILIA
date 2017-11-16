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

function VDJdata = conformGeneGroup(VDJdata, Map, DB)
%Standadize annotations
GrpNum = cell2mat(VDJdata(:, Map.GrpNum));
UnqGrpNum = unique(GrpNum);
for y = 1:length(UnqGrpNum)
    IdxLoc = UnqGrpNum(y) == GrpNum;
    if sum(IdxLoc) <= 1; continue; end %Nothing to correct, single entry.
    Tdata = VDJdata(IdxLoc, :);    
    
    if Map.hSeq > 0 %Heavy chain
        %Standardize the other fields as well (1st one is the root)
        Tdata(:, [Map.hLength(:); Map.hGeneName(:); Map.hGeneNum(:); Map.hDel(:)]) = repmat(Tdata(1, [Map.hLength(:); Map.hGeneName(:); Map.hGeneNum(:); Map.hDel(:)]), size(Tdata, 1), 1);
    end
    
    if Map.lSeq > 0 %Light chain
        %Standardize the other fields as well (1st one is the root)
        Tdata(:, [Map.lLength(:); Map.lGeneName(:); Map.lGeneNum(:); Map.lDel(:)]) = repmat(Tdata(1, [Map.lLength(:); Map.lGeneName(:); Map.lGeneNum(:); Map.lDel(:)]), size(Tdata, 1), 1);
    end
    
    VDJdata(IdxLoc, :) = Tdata;  
end

%Update all data
VDJdata = buildRefSeq(VDJdata, Map, DB, Map.Chain, 'germline', 'first'); %must do first seq of each cluster
VDJdata = updateVDJdata(VDJdata, Map, DB);
