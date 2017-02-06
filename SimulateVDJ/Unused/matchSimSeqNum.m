%matchSimSeqNum will take the simulated data set, and then get the
%annotated dataset sequence order to match with that of the simulated data
%set.


[RefVDJdata,VDJheader,~,~] = openSeqData;
[ShmVDJdata,VDJheader,~,~] = openSeqData;

Seq1 = cell2mat(RefVDJdata(:,H.SeqNumLoc));
Seq2 = cell2mat(ShmVDJdata(:,H.SeqNumLoc));
SeqMap = zeros(size(Seq1,1),1)
for j = 1:size(Seq2,1)
    SeqMap(j) = find(Seq2(j) == Seq1);
end

ShmVDJdata(SeqMap,:) = ShmVDJdata

[FileName,FilePath] = uiputfile('*.xlsx');
xlswrite([FilePath FileName],[VDJheader;ShmVDJdata])
