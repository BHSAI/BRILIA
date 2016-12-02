%copyEquivMatch will use the simulated, singletone unmutated sequences
%equiv matches, and then propogate these equivalent match genes to the
%clonally expanded sequences. To use this:
%1) Find File1 from "generateVDJlib" processed with "findEquivMatch". F
%2) Find File2 from "generateSHMlib"
%3) Run this script copyEquivMatch, looking for the File1 and File2.

[RefVDJdata,NewHeader,~,~] = openSeqData;
[ShmVDJdata,NewHeader,~,~] = openSeqData;

getHeaderVar;

GrpNum = cell2mat(ShmVDJdata(:,GrpNumLoc));
UnqGrpNum = unique(GrpNum);


for y = 1:length(UnqGrpNum)
    IdxLoc = find(UnqGrpNum(y) == GrpNum);
    
    RefIdx = IdxLoc(1);
    RefSeq =ShmVDJdata{RefIdx,RefSeqLoc};
    
    for k = 1:size(RefVDJdata)
        if min(RefSeq == RefVDJdata{k,SeqLoc}) == 1;
            ShmVDJdata(IdxLoc,FamNumLoc) = repmat(RefVDJdata(k,FamNumLoc),length(IdxLoc),1);
            ShmVDJdata(IdxLoc,FamLoc) = repmat(RefVDJdata(k,FamLoc),length(IdxLoc),1);
        end
    end
end

[FileName,FilePath] = uiputfile('*.xlsx');
xlswrite([FilePath FileName],[NewHeader;ShmVDJdata])