%openFastaq will open a fastaq file format, and then return a structure of
%sequence name and format, to be convertable to fasta. This can also be
%saved as an excel or csv file format. Will automatically save a fasta file
%.

% @ERR849860.106 I0VJBII02G7ONM/3
% GGAGACATTGGAGGACTGACAGGAGAAGAGAGACCCA
% +
% GFDGEGEDD11244:<<<<A:99<:88:::::>111>

function SeqData = openFastq

[FileNames, FilePath] = uigetfile('*.fastq','Selection the FASTAQ file','multiselect','on');

if ischar(FileNames)
    FileNames = {FileNames};
end

for j = 1:length(FileNames)
    FileName = FileNames{j};
    FID = fopen([FilePath FileName],'r');
    
    %Determine how many sequences there are first
    NumSeq = 0;
    fseek(FID,0,'eof');
    fmax = ftell(FID);
    frewind(FID);
    
    %LineCount goes 1,2,3,4,1,2,3,4
    LineCount = 1;    
    while feof(FID) == 0
        TextLine = fgetl(FID);
        if LineCount == 1
            NumSeq = NumSeq + 1;            
        end
        LineCount = LineCount + 1;
        if LineCount > 4
            LineCount = 1;
        end
        
        if mod(NumSeq,5000) == 0
            [ftell(FID)/fmax];
        end
    end
    
    SeqData(1:NumSeq) = struct('Header',[],'Sequence',[]);
    
    frewind(FID);
    j = 1;
    LineCount = 1;
    while feof(FID) == 0
        TextLine = fgetl(FID);
        if LineCount == 1
            SeqData(j).Header = TextLine(2:end);
            SeqData(j).Sequence = fgetl(FID);
            LineCount = LineCount + 2;
            j = j+1;
        else
            LineCount = LineCount + 1;
        end
        if LineCount > 4
            LineCount = 1;
        end
        if mod(j,1000) == 0
            [ftell(FID)/fmax];
        end
    end
    
    %Save this as FASTA file for now
    DotLoc = find(FileName == '.');
    SaveName = [FileName(1:DotLoc(end)-1) '.fa'];
    fastawrite([FilePath SaveName],SeqData);   
end


