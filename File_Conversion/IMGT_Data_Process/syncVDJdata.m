%syncVDJdata will take VDJdata1 and VDJdata2, and synce the column that is
%specified according to sequence number match. Mainly used since IMGT does
%not store the other information.
%
%
%  syncVDJdata(SourceVDJdata,DestVDJdata,VDJheader,SyncColNum) will ask first for the orig
%  file, and then the destination file.


function VDJdata2 = syncVDJdata(varargin)
if isempty(varargin)
    [VDJdata1,VDJheader,FileName,FilePath] = openSeqData;
    [VDJdata2,VDJheader,FileName,FilePath] = openSeqData;
else
    VDJdata1 = varargin{1};
    VDJdata2 = varargin{2};
    VDJheader = varargin{3};
%    SyncColNum = varargin{4};
end
H = getHeaderVar(VDJheader);

SeqNum1 = cell2mat(VDJdata1(:,H.SeqNumLoc));
SeqNum2 = cell2mat(VDJdata2(:,H.SeqNumLoc));
[SeqNum12, SeqIdx1, SeqIdx2] = intersect(SeqNum1,SeqNum2); 

VDJdata2(SeqIdx2,H.TemplateLoc) = VDJdata1(SeqIdx1,H.TemplateLoc);

DotLoc = find(FileName == '.');
FileNamePre = [FileName(1:DotLoc(end)-1) 'b'];

%Before saving to xlsx, convert columns with matrix values into char
VDJdata2 = reformatAlignment(VDJdata2,1);
for q = 1:size(VDJdata2,1)
    for w = 1:3
        VDJdata2{q,H.FamNumLoc(w)} = mat2str(VDJdata2{q,H.FamNumLoc(w)});
    end
end

if isempty(varargin)
    %Save to excel or csv file, depending on OS
    if ispc
        xlswrite([FilePath FileNamePre '.xlsx'],cat(1,VDJheader,VDJdata2));
    else
        writeDlmFile(cat(1,VDJheader,VDJdata2),[FilePath FileNamePre '.csv'],'\t');
    end   
end
