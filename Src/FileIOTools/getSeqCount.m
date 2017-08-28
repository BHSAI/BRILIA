%Returns the number of entries in a sequence files
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
    [FID, MSG] = fopen(FullFileName, 'r');
    if FID < 0 
        error('%s: Error opening fasta file [ %s ].\n  %s', mfilename, FullFileName, MSG);
    end
    A = fgetl(FID);
    SeqCount = 0;
    while ischar(A)
        if ~isempty(A) && A(1) == '>'
            SeqCount = SeqCount + 1;
        end
        A = fgetl(FID);
    end
elseif strcmpi(InFileExt, '.fastq')
    [FID, MSG] = fopen(FullFileName, 'r');
    if FID < 0 
        error('%s: Error opening fastq file [ %s ].\n  %s', mfilename, FullFileName, MSG);
    end
    A = fgetl(FID);
    SeqCount = 0;
    while ischar(A)
        if ~isempty(A) && A(1) == '@'
            SeqCount = SeqCount + 1;
            for j = 1:3
                A = fgetl(FID);
                if ~ischar(A)
                    break
                end
            end
        end
        A = fgetl(FID);
    end
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
