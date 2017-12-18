%Returns the number of entries in a sequence files.
%
%  SeqCount = getSeqCount(FullFileName)
%
%  INPUT
%    FullFileName: file name or path + file name of the sequence file
%
%  OUTPUT
%    SeqCount: the number of sequences stored in the sequence file

function SeqCount = getSeqCount(varargin)
%Determine the file type here
if isempty(varargin)
    [InFileName, InFilePath] = uigetfile('*.fa*;*.*sv', 'Select the input sequence file', 'MultiSelect', 'off');
    if isempty(InFileName)
        error('getSeqCount: No file was selected.');
    end
    FullFileName = [InFilePath InFileName];
else
    FullFileName = varargin{1};
    if ~exist(FullFileName, 'file')
        error('%s: Could not fine file [ %s ].', mfilename, FullFileName);
    end
end
[~, ~, InFileExt] = parseFileName(FullFileName);    

%Check to see if there is a file type override
if strcmpi(InFileExt, '.fa') || strcmpi(InFileExt, '.fasta')
    Info = fastainfo(FullFileName);
    SeqCount = Info.NumberOfEntries;
elseif strcmpi(InFileExt, '.fastq')
    Info = fastqinfo(FullFileName);
    SeqCount = Info.NumberOfEntries;
elseif ismember(lower(InFileExt), {'.csv', '.tsv', '.ssv'})
    [FID, MSG] = fopen(FullFileName, 'r');
    if FID < 0
        error('%s: Error opening delimited file [ %s ].\n  %s', mfilename, FullFileName, MSG);
    end
    A = fgetl(FID);
    SeqCount = -1; %For the header, start at -1
    while ischar(A)
        SeqCount = SeqCount + 1;
        A = fgetl(FID);
    end
else
    error('FileType cannot be determined. Make sure file extension is correct');
end
