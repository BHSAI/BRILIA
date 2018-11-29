%convertVDJdata2Fasta will ask  users to select multiple files of the
%VDJdata format, and then convert them to corresponding FASTA files for
%processing. By default, sequence will be named according to the "SeqNum"
%column if  it exists. if not, it will do a simple numerical ordering.

function convertVDJdata2FastaForBaseline(varargin)
[FileNames, FilePath] = uigetfile('*.csv','Select the sequence data files','Multiselect','on');
if ischar(FileNames)
    FileNames = {FileNames};
end

for j = 1:length(FileNames)
    FileName = FileNames{j};
    [VDJdata, VDJheader] = openSeqData([FilePath FileName]);
    DotLoc = find(FileName == '.');
    SaveName = [FilePath, FileName(1:DotLoc(end)-1) 'bsln.fa'];
    saveSeqDataForBaseline(SaveName, VDJdata, VDJheader);
end
