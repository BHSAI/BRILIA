%function findCDR3SpecMut will look for specific mutations that the user
%defined, and return the index location of where that CDR3 mutation has
%occured. Note: Requires that the reference NT be resolved first.

function [MutIdx, MutLoc] = findCDR3SpecMut(varargin)
if length(varargin) == 2
    VDJdata = varargin{1};
    NewHeader = varargin{2};
else
    [VDJdata, NewHeader, ~, ~] = openSeqData([],'eval');
end
AA1 = 'A';
AA2 = 'STV';

CDR3Loc = findHeader(NewHeader,'aminoAcid');
RefLoc = findHeader(NewHeader,'RefSeq');
ClassLoc = findHeader(NewHeader,'Classifier');

MutIdx = zeros(size(VDJdata,1),1);
MutLoc = cell(size(VDJdata,1),1);

for j = 1:size(VDJdata,1)
    [CDR3, ~] = findCDR3(VDJdata(j,:),NewHeader,'RefSeq');
    CDR3Ref = CDR3{1};
    
    CDR3Sam = VDJdata{j,CDR3Loc};

    if isempty(CDR3Sam) || isempty(CDR3Ref)
        disp('No CDR3 Detected for seq');
        continue
    end
    
    %Check for a * in the RefSeq.
    StarLoc = regexp(CDR3Ref,'\*');
    if ~isempty(StarLoc)
        CDR3SamCons = seqconsensus(CDR3Sam);
        for k = 1:length(StarLoc)
            CDR3Ref(StarLoc(k)) = CDR3SamCons(StarLoc(k));
        end
    end
    [~, Ploc] = calcPdist(CDR3Ref,CDR3Sam);
    
    if sum(Ploc) > 0 %Mutation detected
    
        AAref = CDR3Ref(Ploc);
        AAsam = CDR3Sam(Ploc);

        %Decide of this is the mutation you are looking for
        MatchStr1 = regexpi(AAref,['[' AA1 ']*']);
        MatchStr2 = regexpi(AAsam,['[' AA2 ']*']);
    
        if ~isempty(MatchStr1) && ~isempty(MatchStr2)
            MutIdx(j,1) = j;
            MutLoc{j,1} = find(Ploc == 1);
        end
    end    
end

DelLoc = MutIdx == 0;
MutIdx(DelLoc) = [];
MutLoc(DelLoc) = [];


Tdata = VDJdata(MutIdx,:);
CDR3s = Tdata(:,4);
for j = 1:size(Tdata,1)
    MutLocIdx = MutLoc{j};
    CDR3aa = Tdata{j,4};
    CDR3aa(MutLocIdx) = lower(CDR3aa(MutLocIdx));
    CDR3s{j} = CDR3aa;
end
Tdata(:,4) = CDR3s;

%Before saving to xlsx, convert columns with matrix values into char for saving
AddHeader3 = {'vMapNum' 'dMapNum' 'jMapNum'};
SetLoc3 = findHeader(NewHeader,AddHeader3);
for d1 = 1:size(Tdata,1)
    for d2 = 1:3
        Tdata{d1,SetLoc3(d2)} = mat2str(Tdata{d1,SetLoc3(d2)});
    end
end

%Restructure the alignment data
Tdata = reformatAlignment(Tdata,1);
[OutputFilePre, SavePath] = uiputfile('*.*','save as');
if ispc
    xlswrite([OutputFilePre '.cdr3' AA1 '2' AA2 'Filt.xlsx'],cat(1,NewHeader,Tdata));
else
    writeDlmFile(cat(1,NewHeader',Tdata),[OutputFilePre '.Proc.csv'],'\t');
end