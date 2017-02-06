%copyEquivMatch will use the simulated, singletone unmutated sequences
%equiv matches, and then propogate these equivalent match genes to the
%clonally expanded sequences. To use this:
%1) Find File1 from "generateVDJlib" processed with "findEquivMatch". F
%2) Find File2 from "generateSHMlib"
%3) Run this script copyEquivMatch, looking for the File1 and File2.

[RefVDJdata,VDJheader,~,~] = openSeqData;
[ShmVDJdata,VDJheader,~,~] = openSeqData;

H = getHeaderVar(VDJheader);

GrpNum = cell2mat(ShmVDJdata(:,H.GrpNumLoc));
UnqGrpNum = unique(GrpNum);


for y = 1:length(UnqGrpNum)
    IdxLoc = find(UnqGrpNum(y) == GrpNum);
    
    RefIdx = IdxLoc(1);
    RefSeq =ShmVDJdata{RefIdx,H.RefSeqLoc};
    
    for k = 1:size(RefVDJdata)
        if min(RefSeq == RefVDJdata{k,H.SeqLoc}) == 1;
            ShmVDJdata(IdxLoc,H.FamNumLoc) = repmat(RefVDJdata(k,H.FamNumLoc),length(IdxLoc),1);
            ShmVDJdata(IdxLoc,H.FamLoc) = repmat(RefVDJdata(k,H.FamLoc),length(IdxLoc),1);
        end
    end
end

[FileName,FilePath] = uiputfile('*.xlsx');
xlswrite([FilePath FileName],[VDJheader;ShmVDJdata])
