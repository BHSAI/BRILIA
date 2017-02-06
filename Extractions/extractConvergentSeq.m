%extractConvergentSeq will look at multiple sequence files, cluster the
%combined seq using CDR3 info, and then identify unique clusters with
%convergent seq. 

[FileNames,FilePath] = uigetfile('*.xlsx;*.csv','Open the files');
if ischar(FileNames)
    FileNames = {FileNames};
end

%Create a master indext of SeqNum, Seq, CDR3
FileData = cell(length(FileNames),2);
UnqSeqNum = 0; %Going to temporarily assign a unique number to each sequence, during global clustering.

for f = 1:length(FileNames)
    FileData(f,1) = FileNames(f);
    [VDJdata,VDJheader] = openSeqData([FilePath FileNames{f}]);
    H = getHeaderVar(VDJheader);
    
    %Extract the relevent data
    RevData = VDJdata(:,[H.SeqNumLoc H.SeqLoc H.CDR3Loc(1)]);
    UnqSeq = [1:size(RevData,1)]' + UnqSeqNum;
    UnqSeqNum = max(UnqSeq);
    
    FileData = [num2cell(UnqSeq) RevData];
end
