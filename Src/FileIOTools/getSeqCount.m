%Returns the number of entries in a sequence files
function SeqCount = getSeqCount(FullFileName)
%Determine the file type here
if isempty(FullFileName)
    [InFileName, InFilePath] = uigetfile('*.fa*;*.csv', 'Select the input sequence file', 'MultiSelect', 'off');
    if isempty(InFileName);
        error('getSeqCount: No file was selected.');
    end
    FullFileName = [InFilePath InFileName];
end
[~, ~, InFileExt] = parseFileName(FullFileName);    

%Check to see if there is a file type override
if strcmpi(InFileExt, '.fa') || strcmpi(InFileExt, '.fasta')
    Info = fastainfo(FullFileName);
    SeqCount = Info.NumberOfEntries;
elseif strcmpi(InFileExt, '.fastq')
    Info = fastqinfo(FullFileName);
    SeqCount = Info.NumberOfEntries;
elseif strcmpi(InFileExt, '.csv') || strcmpi(InFileExt, '.tsv')
    FID = fopen(FullFileName, 'r');
    A = fgetl(FID);
    SeqCount = -1; %For the header, start at -1
    while ischar(A)
        SeqCount = SeqCount + 1;
        A = fgetl(FID);
    end
else
    error('FileType cannot be determined. Make sure file ext is correct');
end
