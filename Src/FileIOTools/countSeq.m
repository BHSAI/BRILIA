%countSeq returns the number of sequence entries in a sequence files.
%
%  SeqCount = getSeqCount(FileName)
%
%  INPUT
%    FileName: file name of the sequence file (.csv, .fa*). If empty, will
%      ask user to select one.
%
%  OUTPUT
%    SeqCount: the number of sequences stored in the sequence file

function SeqCount = countSeq(FileName)
if nargin == 0 || isempty(FileName)
    [InFileName, InFilePath] = uigetfile('*.fa*;*.*sv', 'Select the input sequence file', 'MultiSelect', 'off');
    assert(ischar(InFileName), '%s: No file was selected.', mfilename);
    FileName = fullfile(InFilePath, InFileName);
else
    assert(exist(FileName, 'file')>0, '%s: Could not find file "%s".', mfilename, FileName);
end

if endsWith(FileName, {'.fa', 'fasta'}, 'ignorecase', true)
    Info = fastainfo(FileName);
    SeqCount = Info.NumberOfEntries;
elseif endsWith(FileName, '.fastq', 'ignorecase', true)
    Info = fastqinfo(FileName);
    SeqCount = Info.NumberOfEntries;
elseif endsWith(FileName, {'.csv', '.tsv', '.ssv'}, 'ignorecase', true)
    [FID, MSG] = fopen(FileName, 'r');
    assert(FID > 0, '%s: Error opening delimited file "%s".\n  %s', mfilename, FileName, MSG);
    TXT = textscan(FID, '%s', 'delimiter', '\n');
    fclose(FID);
    SeqCount = size(TXT{1}, 1) - 1; %Subtract 1 to skip header count.
else
    error('%s: Unrecognized file extension "%s".', mfilename, FileName(find(FileName == '.', 1, 'last'):end));
end

if isempty(SeqCount)
    SeqCount = 0;
end