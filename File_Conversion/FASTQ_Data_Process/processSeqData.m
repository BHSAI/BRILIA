%processSeqData will go through the NGS data from Andrew M. Collins1, Yan
%Wang1, Krishna M. Roskin2, Christopher P. Marquis1 and Katherine J. L.
%Jackson1, PRJEB8745 (www.ebi.ac.uk/ena). Filters out too short sequences,
%and takes out the sequences with a known J. Then cuts and trims them to
%125 bps.

function processSeqData

[FileNames, FilePath] = uigetfile('*.fa','Selection the FASTA file','multiselect','on');

if ischar(FileNames)
    FileNames = {FileNames};
end

[Vmap,Dmap,Jmap] = getCurrentDatabase;

%Obtain the VDJdata header info for output format
[~, ~, StandardData] = xlsread('Headers.xlsx');
NewHeaderLoc = findHeader(StandardData(1,:),'VDJdata');
NewHeader = StandardData(2:end,NewHeaderLoc);
for j = 1:length(NewHeader)
    if isnan(NewHeader{j}); break; end
end
NewHeader(j:end) = [];
NewHeader = NewHeader';

getHeaderVar;
MinLen = 125;

for f = 1:length(FileNames)
    FileName = FileNames{f};
    SeqData = fastaread([FilePath FileName]);

    KeepThis = zeros(size(SeqData))>1;
    for j = 1:length(KeepThis)
        if length(SeqData(j).Sequence) >= MinLen
            SeqF = SeqData(j).Sequence;
            JmatchF = findGeneMatch(SeqF,Jmap,'D'); %We set to 'D' mode since we can have a larger constant region too. Need mid matching.
            SeqR = seqrcomplement(SeqF);
            JmatchR = findGeneMatch(SeqR,Jmap,'D');
            if JmatchF{5}(2) > JmatchR{5}(2)
                Jmatch = JmatchF;
                Seq = SeqF;
            elseif JmatchR{5}(2) > JmatchF{5}(2)
                Jmatch = JmatchR;
                Seq = SeqR;
            else
                if JmatchF{4}(1) > JmatchR{4}(1)
                    Jmatch = JmatchF;
                    Seq = SeqF;
                else
                    Jmatch = JmatchR;
                    Seq = SeqR;
                end                
            end

            JmapNum = Jmatch{1}(1);
            WlocInc = Jmap{JmapNum,end}+2; %Want end of the codon, including + 1. 
            JrefDel = Jmatch{3}(1);
            LeftSeq = Jmatch{4}(1);

            %See if you have enough length
            ValidLen = LeftSeq - JrefDel + WlocInc;
            if  ValidLen >= MinLen && ValidLen <= length(Seq)
                SeqData(j).Sequence = Seq(ValidLen-MinLen+1:ValidLen);
                KeepThis(j) = 1;
            end
        end
        if mod(j,100) == 0;
            [j/length(KeepThis)]
        end
    end
    
    SeqData = SeqData(KeepThis);
    
    DotLoc = find(FileName == '.');
    SaveName = [FileName(1:DotLoc(end)-1) '.Trimmed.fa'];
    SavePath = [cd '\'];
    fastawrite([SavePath SaveName],SeqData);    
    
    
    %Now restructure the data to fit in to excel proper format, using only
    %unique sequences.
    A = struct2cell(SeqData);
    A = squeeze(A)';
    [~, UnqIdx, ~] = unique(A(:,2));
    B = A(UnqIdx,:);
    
    for j = 1:size(B,1)
        SeqOldNum = B{j,1};
        SpaceLoc = regexp(SeqOldNum,'\s','once');
        SeqNum = SeqOldNum(11:SpaceLoc-1);
%        SeqNumCell = regexp(SeqOldNum,'\.(\d*)\s','tokens');
%        B{j,1} = eval(SeqNumCell{1}{1});
        B{j,1} = eval(SeqNum);
    end
    
    SaveName = [FileName(1:DotLoc(end)-1) '.Trimmed.UnqSeq.xlsx'];
    Header = {'SeqNum','nucleotide'};
    xlswrite([SavePath SaveName],[Header;B])
    
    %Go further and ensure the end has TGG ending.
    DelThis = zeros(size(B,1),1) > 1;
    for j = 1:size(B,1)
        SeqEnd = B{j,2}(end-2:end);
        if ~strcmpi(SeqEnd,'TGG')
           DelThis(j) = 1;
        end
    end
    B(DelThis,:) = [];
    
    %Now look productive only
    DelThis = zeros(size(B,1),1) > 1;
    for j = 1:size(B,1)
        AAseq = nt2aa(B{j,2},'frame',3,'ACGTonly','false');
        if ~isempty(strfind(AAseq,'*'))
           DelThis(j) = 1;
        end
    end
    B(DelThis,:) = [];
    
    SaveName = [FileName(1:DotLoc(end)-1) '.Trimmed.UnqSeq.Prod.xlsx'];
    Header = {'SeqNum','nucleotide'};
    xlswrite([SavePath SaveName],[Header;B])    
    
    %Sequences are sorted for now. Look for pairwise distance to look for
    %sudden increases in mismatche. use these to segment the file, around
    %the 500 seq mark.
    DiffMat = zeros(size(B,1),1);
    for j = 2:size(DiffMat,1)
        DiffMat(j) = sum(B{j-1,2} ~= B{j,2});
    end
    SegLocT = find(DiffMat>20);
    SegSize = 1000;
    SegLoc = zeros(ceil(size(B,1)/SegSize),1);
    j = 1;
    while SegSize*j <= size(B,1)
        SeqLocTT = find(SegLocT > SegSize*j);
        SegLoc(j) = SegLocT(SeqLocTT(1));
        j = j+1;
    end
    SegLoc(end) = size(B,1);    
    
    S1 = 1;
    for f = 1:length(SegLoc)
        S2 = SegLoc(f)-1; %Remember, segments define start of sequence.
        SaveName = sprintf('%s_%0.0d-%0.0d.xlsx',FileName(1:DotLoc(end)-1),S1,S2);
        Header = {'SeqNum','nucleotide'};
        xlswrite([SavePath SaveName],[Header;B(S1:S2,:)])            
        S1 = S2+1;
    end
end