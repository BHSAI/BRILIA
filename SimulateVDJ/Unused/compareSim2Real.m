%compareSim2Real will select the File1-Truth and File2-Annotated data, and
%return the error and accuracy statistics.

function compareSim2Real(varargin)
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
SeqLoc1 = H.SeqLoc;
SeqNumLoc1 = H.SeqNumLoc;
GrpNumLoc1 = H.GrpNumLoc;
RefSeqLoc1 = H.RefSeqLoc;
FamLoc1 = H.FamLoc;
FamNumLoc1 = H.FamNumLoc;

[VDJdata2, VDJheader, FileName2, FilePath2] = openSeqData(FullPath2);
H = getHeaderVar(VDJheader);
SeqLoc2 = H.SeqLoc;
SeqNumLoc2 = H.SeqNumLoc;
GrpNumLoc2 = H.GrpNumLoc;
RefSeqLoc2 = H.RefSeqLoc;
FamLoc2 = H.FamLoc;
FamNumLoc2 = H.FamNumLoc;

VDJdata1 = removeNAN(VDJdata1);
VDJdata2 = removeNAN(VDJdata2);


%For the 2nd dataset, you want to ensure all REFSEQ is set to the ROOT's
%REFSEQ. THat is because BRILIA returns parent-child seq, when we want to
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




%Selecting first entry only
for j = 1:size(VDJdata2,1)
    if ~isempty(VDJdata2{j,FamNumLoc2(1)});
        VDJdata2{j,FamNumLoc2(1)} = VDJdata2{j,FamNumLoc2(1)}(1); 
    end
    if ~isempty(VDJdata2{j,FamNumLoc2(2)});
        VDJdata2{j,FamNumLoc2(2)} = VDJdata2{j,FamNumLoc2(2)}(1); 
    end
    if ~isempty(VDJdata2{j,FamNumLoc2(3)});
        VDJdata2{j,FamNumLoc2(3)} = VDJdata2{j,FamNumLoc2(3)}(1); 
    end
end

%Extract the unique sequence numbers
SeqNums1 = cell2mat(VDJdata1(:,SeqNumLoc1));
SeqNums2 = cell2mat(VDJdata2(:,SeqNumLoc2));

%Pick out the same sequences
[~, Idx2, Idx1] = intersect(SeqNums2,SeqNums1);
VDJdata1 = VDJdata1(Idx1,:);
VDJdata2 = VDJdata2(Idx2,:);

%Compare side by side, REF SEQ to REAL SEQ
MisMatchCt = zeros(size(VDJdata1,1),1);
for j = 1:size(VDJdata1,1)
    RefSeq1 = VDJdata1{j,RefSeqLoc1};
    RefSeq2 = VDJdata2{j,RefSeqLoc2};
    if length(RefSeq1) ~= length(RefSeq2)
        MisMatchCt(j) = -1;
        continue
    end
    if isempty(RefSeq1) || isempty(RefSeq2)
        MisMatchCt(j) = -1;
        continue
    end
    if VDJdata1{j,H.FunctLoc} == 'N' || VDJdata2{j,H.FunctLoc} == 'N'
        MisMatchCt(j) = -1;
        continue
    end
    MisMatchCt(j) = sum(RefSeq1~=RefSeq2);
end
MisMatchCt(MisMatchCt<0) = [];
Edges = [0:max(MisMatchCt)]';
N = histc(MisMatchCt,Edges);

%Determine gene name match score if we allow equivalent matches
VgeneStat = zeros(size(VDJdata1,1),1); 
DgeneStat = zeros(size(VDJdata1,1),1);
JgeneStat = zeros(size(VDJdata1,1),1);

for j = 1:size(VDJdata1,1)   
    %Determine how often the gene number is matched correctly
    Vname1 = VDJdata1{j,FamLoc1(1)};
    Dname1 = VDJdata1{j,FamLoc1(2)};
    Jname1 = VDJdata1{j,FamLoc1(3)};
    Vname2 = VDJdata2{j,FamLoc2(1)};
    Dname2 = VDJdata2{j,FamLoc2(2)};
    Jname2 = VDJdata2{j,FamLoc2(3)};
    
    %Remove the "0" if it is located after IGHX
    V1zero = regexp(Vname1,'IGH[VDJ]')+4;
    V1del = V1zero(Vname1(V1zero) == '0');
    Vname1(V1del) = [];   
    D1zero = regexp(Dname1,'IGH[VDJ]')+4;
    D1del = D1zero(Dname1(D1zero) == '0');
    Dname1(D1del) = [];        
    J1zero = regexp(Jname1,'IGH[VDJ]')+4;
    J1del = J1zero(Jname1(J1zero) == '0');
    Jname1(J1del) = [];
    
    V2zero = regexp(Vname2,'IGH[VDJ]')+4;
    V2del = V2zero(Vname2(V2zero) == '0');
    Vname2(V2del) = [];   
    D2zero = regexp(Dname2,'IGH[VDJ]')+4;
    D2del = D2zero(Dname2(D2zero) == '0');
    Dname2(D2del) = [];        
    J2zero = regexp(Jname2,'IGH[VDJ]')+4;
    J2del = J2zero(Jname2(J2zero) == '0');
    Jname2(J2del) = [];
    
    %Replace *'s with letter, since this causes issues with regexp.
    Vname1 = strrep(Vname1,'*','u');
    Vname1 = regexp(Vname1,'\|','split');
    Vname2 = strrep(Vname2,'*','u');
    ColonLoc = find(Vname2 == ':');
    if ~isempty(ColonLoc)
        Vname2(1:ColonLoc+1) = [];
    end
    for q = 1:length(Vname1)
        if sum(strfind(Vname2,Vname1{q})) > 0
            VgeneStat(j) = 1;
            break
        end
    end
    
    %Replace *'s with letter, since this causes issues with regexp.
    Dname1 = strrep(Dname1,'*','u');
    Dname1 = regexp(Dname1,'\|','split');
    Dname2 = strrep(Dname2,'*','u');
    ColonLoc = find(Dname2 == ':');
    if ~isempty(ColonLoc)
        Dname2(1:ColonLoc+1) = [];
    end
    for q = 1:length(Dname1)
        if sum(strfind(Dname2,Dname1{q})) > 0
            DgeneStat(j) = 1;
            break
        end
    end
    
    %Replace *'s with letter, since this causes issues with regexp.
    Jname1 = strrep(Jname1,'*','u');
    Jname1 = regexp(Jname1,'\|','split');
    Jname2 = strrep(Jname2,'*','u');
    ColonLoc = find(Jname2 == ':');
    if ~isempty(ColonLoc)
        Jname2(1:ColonLoc+1) = [];
    end
    for q = 1:length(Jname1)
        if sum(strfind(Jname2,Jname1{q})) > 0
            JgeneStat(j) = 1;
            break
        end
    end
end

%Group identity match by unique D gene of ground truth
Vacc = sum(VgeneStat);%;/size(VDJdata1,1);
Dacc = sum(DgeneStat);%;/size(VDJdata1,1);
Jacc = sum(JgeneStat);%;/size(VDJdata1,1);
VDJacc = sum(VgeneStat & DgeneStat & JgeneStat);

%Setup the output file

%Determine if the initial database is for equiv, or absolute match
if sum(strfind(FileName1,'Equiv')) > 0
    SheetPre = 'Equiv_'; 
else
    SheetPre = 'Abs_';
end

CellHeader1 = {'MisMatch_Count' 'Frequency'};
CellData1 = cell(length(N),2);
CellData1(1:end,1) = num2cell(Edges);
CellData1(1:end,2) = num2cell(N);

CellHeader2 = {'Gene_Segment' 'Match_Accuracy'};
CellData2 = cell(4,2);
CellData2(1:end,1) = {'Vgene'; 'Dgene'; 'Jgene'; 'VDJgene'};
CellData2(1:end,2) = num2cell([Vacc; Dacc; Jacc; VDJacc]);

xlswrite([FilePath2 FileName2],[CellHeader1; CellData1],[SheetPre 'MisMatchCt']);
xlswrite([FilePath2 FileName2],[CellHeader2; CellData2],[SheetPre 'GeneFamMatch']);
