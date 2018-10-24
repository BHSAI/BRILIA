%convertVDJdata2Fasta will ask  users to select multiple files of the
%VDJdata format, and then convert them to corresponding FASTA files for
%processing. By default, sequence will be named according to the "SeqNum"
%column if  it exists. if not, it will do a simple numerical ordering.

function convertVDJdata2Fasta(varargin)
[FileNames, FilePath] = uigetfile('*.*sv','Select the BRILIA output data files','Multiselect','on');
if ischar(FileNames)
    FileNames = {FileNames};
end

for f = 1:length(FileNames)
    FileName = FileNames{f};
    [VDJdata, VDJheader] = openSeqData(fullfile(FilePath, FileName));    
    [H, L, Chain] = getAllHeaderVar(VDJheader);
    SeqNumLoc = H.SeqNumLoc;
    switch Chain
        case 'H'
            SeqLoc = H.SeqLoc;
        case 'L'
            SeqLoc = L.SeqLoc;
        case 'HL'
            SeqLoc = [H.SeqLoc L.SeqLoc];
        otherwise
            continue;
    end
    
    SeqNum = VDJdata(:, SeqNumLoc);
    for c = 1:length(SeqLoc)
        Seq = VDJdata(:, SeqLoc(c));
        for j = 1:length(SeqNum)
            SeqNum{j} = ['Seq_' num2str(SeqNum{j})];
        end
        DataStruct = cell2struct([SeqNum Seq], {'Header' 'Sequence'}, 2);
        DotLoc = find(FileName == '.');
        SaveName = [FilePath FileName(1:DotLoc(end)-1) '_' Chain(c) '.fa'];
        fastawrite(SaveName, DataStruct)
        clear DataStruct Seq 
    end
    clear VDJdata VDJheader SeqNum 
end
