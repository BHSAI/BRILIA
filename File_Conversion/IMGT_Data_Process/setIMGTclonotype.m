%setIMGTclonotype will group VDJdata entries according to IMGT's definition
%of a clone, which is unique V,D,J genes, same CDR3 sequence. HOWEVER, we
%do allow CDR3 seq to vary by N number of amino acids or nts.
%
%  VDJdata = setIMGTclonotype(VDJdata,VDJheader,Alphabet,HamDist) 
%    Alphabet = 'nt' or 'aa'
%    HamDist = hamming distance, based on Alphabet, of the CDR3 region

function VDJdata = setIMGTclonotype()%Alphabet,HamDist)
HamDist = 0.03;
Alphabet = 'nt';
[VDJdata,VDJheader,FileName,FilePath] = openSeqData;
VDJdata = removeNAN(VDJdata);

DotLoc = find(FileName == '.');
FileNamePre = sprintf('%s_%sHam%d%%',FileName(1:DotLoc-1),Alphabet,HamDist*100);

H = getHeaderVar(VDJheader);

%Removing VDJ's without fully resolved V D and J.
DelThese = zeros(size(VDJdata,1),1) > 1;
for j = 1:size(VDJdata,1)
    if VDJdata{j,H.FamNumLoc(1)} == 0
        DelThese(j) = 1;
        continue
    elseif VDJdata{j,H.FamNumLoc(2)} == 0
        DelThese(j) = 1;
        continue
    elseif VDJdata{j,H.FamNumLoc(3)} == 0
        DelThese(j) = 1;
        continue
    elseif size(VDJdata{j,H.SeqLoc},2) ~= 125
        DelThese(j) = 1;
        continue
    end  
end

%cluster VDJdata the same way
VDJdata(DelThese,:) = [];
VDJdata = clusterGeneIMGT(VDJdata,VDJheader);

% %Just in case the input has multiple solutions
% Vcell = VDJdata(:,H.FamNumLoc(1)); %Temp store V fam num
% Dcell = VDJdata(:,H.FamNumLoc(2)); %Temp store D fam num
% Jcell = VDJdata(:,H.FamNumLoc(3)); %Temp store J fam num
% for j = 1:size(VDJdata)
%     VDJdata{j,H.FamNumLoc(1)} = VDJdata{j,H.FamNumLoc(1)}(1); 
%     VDJdata{j,H.FamNumLoc(2)} = VDJdata{j,H.FamNumLoc(2)}(1); 
%     VDJdata{j,H.FamNumLoc(3)} = VDJdata{j,H.FamNumLoc(3)}(1); 
% end
% 
% %Subdivide the clusters based on VDJ uniqueness
% UnqGrp = cell2mat(VDJdata(:,H.GrpNumLoc));
% [~,~,UnqVDJidx] = unique(cell2mat(VDJdata(:,H.FamNumLoc)),'rows');
% [~,~,CloneNum] = unique([UnqGrp UnqVDJidx],'rows') ;
% [CloneNum, CloneIdx] = sort(CloneNum);
% 
% %Just in case the input has multiple solutions
% VDJdata(:,H.FamNumLoc) = [Vcell Dcell Jcell]; %Replace original numbering before resorting
% VDJdata = VDJdata(CloneIdx,:);
% VDJdata(:,H.GrpNumLoc) = num2cell(CloneNum);
%  
%Before saving to xlsx, convert columns with matrix values into char
VDJdataSave = reformatAlignment(VDJdata,1);
for q = 1:size(VDJdata,1)
    for w = 1:3
        VDJdataSave{q,H.FamNumLoc(w)} = mat2str(VDJdataSave{q,H.FamNumLoc(w)});
    end
end

%Save to excel or csv file, depending on OS
if ispc
    xlswrite([FilePath FileNamePre '.xlsx'],cat(1,VDJheader,VDJdataSave));
else
    writeDlmFile(cat(1,VDJheader,VDJdataSave),[FilePath FileNamePre '.csv'],'\t');
end
