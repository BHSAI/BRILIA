%specialized analysis for Ilja's analysis on germline vs maximum diverged
%sequence comparison.
%  1) Takes the germline for each clone --> save to file1
%  2) Take the maximum hamming distance separation from germline clone -->
%  save to file2
%  3) Converts both files into Adaptive format.

[FileNames, FilePath] = uigetfile('*.xlsx','Selection files','multiselect','on');
if ischar(FileNames)
    FileNames = {FileNames};
end

for f = 1:length(FileNames)
    [VDJdata,VDJheader,FileName,FilePath] = openSeqData([FilePath FileNames{f}]);
    H = getHeaderVar(VDJheader);
    
    %Removing the non-functional sequences, that has stop codon or
    %out-of-frame junctions. Pseudogenes?
    DelThis = char(VDJdata(:,H.FunctLoc)) == 'N';
    VDJdata(DelThis,:) = [];
    
    %Identify clonal groups and extract germSeq, maxHamDistSeq.
    GrpNum = cell2mat(VDJdata(:,H.GrpNumLoc));
    UnqGrpNum = unique(GrpNum);
    GermData = cell(length(UnqGrpNum),size(VDJdata,2)); %Stores the germline sequence
    DevData = cell(length(UnqGrpNum),size(VDJdata,2)); %Stores the maximum hamming distance one, tied broken with template count
    
    for y = 1:length(UnqGrpNum)
        %Find the germline sequence for each clonal group
        IdxLoc = find(UnqGrpNum(y) == GrpNum);
        RootLoc = IdxLoc(1);
        GermData(y,:) = VDJdata(RootLoc,:);
        GermData(y,H.SeqLoc) = GermData(y,H.RefSeqLoc); %need to make sure Seq is the RefSeq;
        
        %Find the maximum distance seq
        RefSeq = GermData{y,H.SeqLoc};
        DistTrack = zeros(length(IdxLoc),3); %[HamDist TempCount IdxLoc]
        for k = 1:length(IdxLoc)
            CurSeq = VDJdata{IdxLoc(k),H.SeqLoc};
            TempCt = VDJdata{IdxLoc(k),H.TemplateLoc};
            HamDist = sum(CurSeq ~= RefSeq);
            DistTrack(k,:) = [HamDist TempCt IdxLoc(k)];
        end
        DistTrack = sortrows(DistTrack,[-1 -2]); % Searching for highest deviation, and larger template count.
        MaxLoc = DistTrack(1,end);
        DevData(y,:) = VDJdata(MaxLoc,:);
    end
    
    %Rebuild SHM mutations and classifiers for the GermData
    GermData = appendMutCt(GermData,VDJheader);
    GermData = makeClassifier(GermData,VDJheader);
    GermData = buildRefSeq(GermData,VDJheader,'single');
    GermData = buildVDJalignment(GermData,VDJheader);
    
    %Select output file name
    DotLoc = find(FileName == '.');
    DotLoc = DotLoc(end);
    if ispc
        FileName1 = [FileName(1:DotLoc) 'Germline.xlsx'];
        FileName2 = [FileName(1:DotLoc) 'MaxDiver.xlsx'];
    else
        FileName1 = [FileName(1:DotLoc) 'Germline.csv'];
        FileName2 = [FileName(1:DotLoc) 'MaxDiver.csv'];
    end

    %Converting matrix to text for excel writing
    for q = 1:size(GermData,1)
        for w = 1:3
            GermData{q,H.FamNumLoc(w)} = mat2str(GermData{q,H.FamNumLoc(w)});
            DevData{q,H.FamNumLoc(w)} = mat2str(DevData{q,H.FamNumLoc(w)});
        end
    end

    %Save to excel or csv file, depending on OS
    if ispc
        xlswrite([FilePath FileName1],cat(1,VDJheader,GermData));
        xlswrite([FilePath FileName2],cat(1,VDJheader,DevData));
    else
        writeDlmFile(cat(1,VDJheader,GermData),[FilePath FileName1],'\t');
        writeDlmFile(cat(1,VDJheader,DevData),[FilePath FileName2],'\t');
    end
    
    %Convert to AdapFile, but remove the header
    convertVDJdata2Adaptive(FileName1,FilePath,'NoHeader');
    convertVDJdata2Adaptive(FileName2,FilePath,'NoHeader');
    
end
 
