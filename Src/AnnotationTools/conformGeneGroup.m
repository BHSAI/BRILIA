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
%
function VDJdata = conformGeneGroup(VDJdata, Map, DB)
if isempty(VDJdata); return; end

if ~iscell(VDJdata{1}) %VDJdata is not spliced
    %Standadize annotations
    GrpNum = cell2mat(VDJdata(:, Map.GrpNum));
    UnqGrpNum = unique(GrpNum);
    for y = 1:length(UnqGrpNum)
        IdxLoc = UnqGrpNum(y) == GrpNum;
        if sum(IdxLoc) <= 1; continue; end %Nothing to correct, single entry.
        Tdata = VDJdata(IdxLoc, :);    
        if Map.hSeq > 0 %Heavy chain. %Standardize the other fields as well (1st one is the root)
            Tdata(:, [Map.hLength(:); Map.hGeneName(:); Map.hGeneNum(:); Map.hDel(:)]) = repmat(Tdata(1, [Map.hLength(:); Map.hGeneName(:); Map.hGeneNum(:); Map.hDel(:)]), size(Tdata, 1), 1);
        end
        if Map.lSeq > 0 %Light chain. %Standardize the other fields as well (1st one is the root)
            Tdata(:, [Map.lLength(:); Map.lGeneName(:); Map.lGeneNum(:); Map.lDel(:)]) = repmat(Tdata(1, [Map.lLength(:); Map.lGeneName(:); Map.lGeneNum(:); Map.lDel(:)]), size(Tdata, 1), 1);
        end
        VDJdata(IdxLoc, :) = Tdata;  
    end
    VDJdata = buildRefSeq(VDJdata, Map, DB, Map.Chain, 'germline', 'first'); %must do first seq of each cluster
    VDJdata = updateVDJdata(VDJdata, Map, DB);
else %VDJdata is spliced for parfor
    parfor y = 1:length(VDJdata)
        if isempty(VDJdata{y}); pause; end
        if Map.hSeq > 0 %#ok<PFBNS> %Heavy chain. %Standardize the other fields as well (1st one is the root)
            VDJdata{y}(:, [Map.hLength(:); Map.hGeneName(:); Map.hGeneNum(:); Map.hDel(:)]) = repmat(VDJdata{y}(1, [Map.hLength(:); Map.hGeneName(:); Map.hGeneNum(:); Map.hDel(:)]), size(VDJdata{y}, 1), 1);
        end
        if Map.lSeq > 0 %Light chain. %Standardize the other fields as well (1st one is the root)
            VDJdata{y}(:, [Map.lLength(:); Map.lGeneName(:); Map.lGeneNum(:); Map.lDel(:)]) = repmat(VDJdata{y}(1, [Map.lLength(:); Map.lGeneName(:); Map.lGeneNum(:); Map.lDel(:)]), size(VDJdata{y}, 1), 1);
        end
        VDJdata{y} = buildRefSeq(VDJdata{y}, Map, DB, Map.Chain, 'germline', 'first'); %must do first seq of each cluster
        VDJdata{y} = updateVDJdata(VDJdata{y}, Map, DB);
    end
end