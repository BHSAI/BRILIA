%filterProdVDJ will take an excel file, translate the Seq, and select only
%sequences that are coding (no stop codons);

function filterNonprodVDJ(varargin)
%Open file.
[VDJdata, NewHeader, FileName, FilePath] = openSeqData([],'false');
AACol = findHeader(NewHeader,'aminoacid');
AAseq = VDJdata(:,AACol);

%Remove any CDR3's with *, or is empty, or lacks the conserved C or F/W. 
ProdSeqLoc = zeros(size(AAseq));
for j = 1:size(AAseq,1)
    if sum(isnan(AAseq{j})) >= 1
        continue;
    end
    if isempty(AAseq{j})
        continue;
    end
    if sum(AAseq{j} == '*') == 0 && AAseq{j}(1,1) == 'C' && AAseq{j}(1,end) == 'W'
        ProdSeqLoc(j) = 1;
    end
end

%Condense data and save
ProdData = VDJdata(ProdSeqLoc==1,:);
NonprodData = VDJdata(ProdSeqLoc==0,:);

%Save the sorted file
DotLoc = find(FileName == '.');
OutputFilePre = FileName(1:DotLoc(end));

%Fix the group number for productive, if it exist.
GrpCol = findHeader(NewHeader,'GroupNum'); %Column where group number is often stored
if ~isempty(ProdData) && GrpCol ~= 0
    if ~isempty(ProdData{1,GrpCol})
        OldGroupNum = cell2mat(ProdData(:,GrpCol));
        [~, ~, NewGrp] = unique(OldGroupNum);
        [NewGrp, ResortIdx] = sort(NewGrp);
        ProdData = ProdData(ResortIdx,:);
        ProdData(:,GrpCol) = num2cell(NewGrp);
    end
    
    if ispc %Write to Excel if it's a pc
        xlswrite([FilePath OutputFilePre 'prod.xlsx'],cat(1,NewHeader,ProdData));
    else %Write to csv, tab-delimited file
        writeDlmFile(cat(1,NewHeader,ProdData),[FilePath OutputFilePre 'prod.csv'],'\t');
    end 
end

%Fix the group number for Nonproductive, if it exist.
if ~isempty(NonprodData) && GrpCol ~= 0
    if ~isempty(NonprodData{1,GrpCol})
        OldGroupNum = cell2mat(NonprodData(:,GrpCol));
        [~, ~, NewGrp] = unique(OldGroupNum);
        [NewGrp, ResortIdx] = sort(NewGrp);
        NonprodData = NonprodData(ResortIdx,:);
        NonprodData(:,GrpCol) = num2cell(NewGrp);
    end
    
    if ispc %Write to Excel if it's a pc
        xlswrite([FilePath OutputFilePre 'nonprod.xlsx'],cat(1,NewHeader,NonprodData));
    else %Write to csv, tab-delimited file
        writeDlmFile(cat(1,NewHeader,NonprodData),[FilePath OutputFilePre 'nonprod.csv'],'\t');
    end
end