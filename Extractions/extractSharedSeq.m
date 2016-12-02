%this script looks for shared CDR3s, and then extracts the groups that have
%them, saves them into one file.
[SeqData,SeqHeader,FileName,FilePath] = openSeqData();%['C:\Users\dlee\Desktop\Matlab Script\Dev9\Specialized_analysis\SharedSeq.xlsx']);
SharedCDR3Loc = findHeader(SeqHeader,'aminoAcid');
SharedCDR3 = SeqData(:,SharedCDR3Loc);

[FileNames, FilePath] = uigetfile('*.xlsx','Select Files to look for shared seq','multiselect','on');
if ~iscell(FileNames)
    FileNames = {FileNames};
end

for f = 1:length(FileNames)
    FileName = FileNames{f};
    
    [VDJdata,NewHeader,~,~] = openSeqData([FilePath FileName]);
    getHeaderVar;

    GrpNum = cell2mat(VDJdata(:,GrpNumLoc));
    GrpNumUnq = unique(GrpNum);
    KeepGrpNum = zeros(length(GrpNumUnq),1) > 1;
    NumOfElem = 0;
    
    %Identifying all group that have same CDR3 as the shared seq
    for j = 1:size(VDJdata,1)
        CDR3 = VDJdata{j,CDR3Loc(1)};
        CurGrpNum = VDJdata{j,GrpNumLoc};
        if KeepGrpNum(CurGrpNum == GrpNumUnq) == 1
            continue
        end
        
        for k = 1:size(SharedCDR3,1)
            CDR3share = SharedCDR3{k};
            if length(CDR3share) == length(CDR3)
                if strcmpi(CDR3share,CDR3)
                    KeepGrpNum(GrpNumUnq == CurGrpNum) = 1;
                    NumOfElem = NumOfElem + sum(GrpNum == CurGrpNum);
                    break
                end
            end
        end
    end
    
    %Extract sequence with the same group and CDR3
    SelectGrpNum = GrpNumUnq(KeepGrpNum);
    NewVDJdata = cell(NumOfElem,length(NewHeader));
    r = 1;
    for q = 1:length(SelectGrpNum)
        IdxLoc = find(GrpNum == SelectGrpNum(q));
        NewVDJdata(r:r+length(IdxLoc)-1,:) = VDJdata(IdxLoc,:);
        r = r+length(IdxLoc);
    end
    %Before saving to xlsx, convert columns with matrix values into char
    NewVDJdata = reformatAlignment(NewVDJdata,1);
    for q = 1:size(NewVDJdata,1)
        for w = 1:3
            NewVDJdata{q,FamNumLoc(w)} = mat2str(NewVDJdata{q,FamNumLoc(w)});
        end
    end
    
    %Saving just the sharedSeq
    SaveName = ['SharedSeq_' FileName];
    xlswrite(SaveName,[NewHeader;NewVDJdata]);
end
        