%compareSim2Real will select the File1-Truth and File2-Annotated data, and
%return the error and accuracy statistics.

function compareSim2Real_SHM(varargin)
P = inputParser;
addOptional(P,'FullPath1',[],@ischar);
addOptional(P,'FullPath2',[],@ischar);
parse(P,varargin{:});
FullPath1 = P.Results.FullPath1;
FullPath2 = P.Results.FullPath2;

%Select the truth file
if isempty(FullPath1)
    [FileName1, FilePath1] = uigetfile('*.xlsx','Select truth file');
    FullPath1 = [FilePath1 FileName1];
end
if isempty(FullPath2)
    [FileName2, FilePath2] = uigetfile('*.xlsx','Select tested annoted file');
    FullPath2 = [FilePath2 FileName2];
end

%Select the annotated file
[VDJdata1, VDJheader, FileName1, FilePath1] = openSeqData(FullPath1);
H = getHeaderVar(VDJheader);
VDJdata1 = removeNAN(VDJdata1);
SeqLoc1 = H.SeqLoc;
SeqNumLoc1 = H.SeqNumLoc;
GrpNumLoc1 = H.GrpNumLoc;
RefSeqLoc1 = H.RefSeqLoc;

[VDJdata2, VDJheader, FileName2, FilePath2] = openSeqData(FullPath2);
H = getHeaderVar(VDJheader);
VDJdata2 = removeNAN(VDJdata2);
SeqLoc2 = H.SeqLoc;
SeqNumLoc2 = H.SeqNumLoc;
GrpNumLoc2 = H.GrpNumLoc;
RefSeqLoc2 = H.RefSeqLoc;

%For the 2nd dataset, you want to ensure all REFSEQ is set to the ROOT's
%REFSEQ. That is because BRILIA returns parent-child seq, when we want to
%compare germline-child relations.
GrpNum1 = cell2mat(VDJdata1(:,GrpNumLoc1));
UnqGrpNum1 = unique(GrpNum1);
for y = 1:length(UnqGrpNum1)
    IdxLoc = find(UnqGrpNum1(y) == GrpNum1);
    VDJdata1(IdxLoc,RefSeqLoc2) = VDJdata1(IdxLoc(1),RefSeqLoc1);
end

GrpNum2 = cell2mat(VDJdata2(:,GrpNumLoc2));
UnqGrpNum2 = unique(GrpNum2);
for y = 1:length(UnqGrpNum2)
    IdxLoc = find(UnqGrpNum2(y) == GrpNum2);
    VDJdata2(IdxLoc,RefSeqLoc2) = VDJdata2(IdxLoc(1),RefSeqLoc2);
end


%Resort true VDJ based on grp, and then seq num. That way, we know which
%has the 0,5,10,15,20,25 SHMs. Going to separate results by these numbers,
%to calculates TP FN FP TN.
VDJdata1 = sortrows(VDJdata1,[GrpNumLoc1 SeqNumLoc1]);

%Iteratively find groups
GrpNum1 = cell2mat(VDJdata1(:,GrpNumLoc1));
SeqNum2 = cell2mat(VDJdata2(:,SeqNumLoc2));
UnqGrpNum1 = unique(GrpNum1);
StatData = zeros(6,4);
for y = 1:length(UnqGrpNum1)
    IdxLoc = find(UnqGrpNum1(y) == GrpNum1);
    CDR3start = VDJdata1{IdxLoc(1),H.CDR3Loc(3)};
    CDR3end = VDJdata1{IdxLoc(1),H.CDR3Loc(4)};
    for k = 1:length(IdxLoc)
        Seq1 = VDJdata1{IdxLoc(k),SeqLoc1};
        RefSeq1 = VDJdata1{IdxLoc(k),RefSeqLoc1};
        SeqNum1 = VDJdata1{IdxLoc(k),SeqNumLoc1};
        SeqIdx2 = find(SeqNum2 == SeqNum1);
        Seq2 = VDJdata2{SeqIdx2,SeqLoc2};
        RefSeq2 = VDJdata2{SeqIdx2,RefSeqLoc2};
        
        if isempty(RefSeq1) || isempty(RefSeq2) || isempty(SeqIdx2) || length(RefSeq1) ~= length(RefSeq2);
            continue
        end
        
        %Focus only on CDR3 length
        Seq1 = Seq1(CDR3start:CDR3end);
        Seq2 = Seq2(CDR3start:CDR3end);
        RefSeq1 = RefSeq1(CDR3start:CDR3end);
        RefSeq2 = RefSeq2(CDR3start:CDR3end);
        
        %Find location of real and cur SHMs
        RealSHM = RefSeq1 ~= Seq1;
        CurrSHM = RefSeq2 ~= Seq2;
        CombSHM = RealSHM + CurrSHM;
        EvalLoc = find(CombSHM > 0);
        
        TP = 0; %Found correct TP
        FP = 0; %Found SHM that did not exist
        for q = 1:length(EvalLoc)
            RealX0 = RefSeq1(EvalLoc(q));
            RealX1 = Seq1(EvalLoc(q));
            CurrX0 = RefSeq2(EvalLoc(q));
            CurrX1 = Seq2(EvalLoc(q));
            
            if strcmpi([RealX0 RealX1],[CurrX0 CurrX1])
                TP = TP + 1;
            end
            
            if strcmpi(RealX0,RealX1) && ~strcmpi(CurrX0,CurrX1)
                FP = FP + 1;
            end
            
        end
        FN = sum(RealSHM) - TP;
        TN = length(RefSeq1) - TP - FP - FN;
        
        StatData(k,:) = StatData(k,:) + [TP TN FP FN];
    end    
end

Summary = [[0:5:25]' sum(StatData,2)/125 StatData];
SumHeader = {'TrueSHMct' 'SeqCoverage' 'TP' 'TN' 'FP' 'FN'};
WriteData = [SumHeader; num2cell(Summary)];
xlswrite([FilePath2 FileName2],WriteData,'SHMstats');
