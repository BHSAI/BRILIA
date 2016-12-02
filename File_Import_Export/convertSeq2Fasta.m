%convertSeq2Fasta will ask  users to select multiple files containing
%either the VDJdata or Adap file format, and then convert them to
%corresponding FASTA files for processing. By default, sequence will be
%named according to the "SeqNum" column if exists. if not, it will do a
%simple numerical ordering.

function convertSeq2Fasta(varargin)

[FileNames, FilePath] = uigetfile('*.xlsx;*.csv','Select the sequence data files','Multiselect','on');

if ischar(FileNames)
    FileNames = {FileNames};
end

for j = 1:length(FileNames)
    FileName = FileNames{j};
    [VDJdata, NewHeader, ~, ~] = openSeqData([FilePath FileName]);    
    getHeaderVar;
    
    %Assign seq header name if none is inputted
    if isempty(varargin)
        if SeqLoc == 0
            error('Error: Could not find the nucleotide column');
        else
            Seq = VDJdata(:,SeqLoc);
        end

        if SeqNumLoc == 0
            disp('Warning: Could not find the SeqNum column - making own numbering');
            SeqNum = num2cell([1:size(VDJdata,1)]');
        else
            SeqNum = VDJdata(:,SeqNumLoc);
        end

        Header = cell(size(VDJdata,1),1);
        for k = 1:size(Header,1)
            Header{k} = sprintf('Seq_%d',SeqNum{k});
        end
    else
        HeaderLoc = varargin{1};
        SeqLoc = varargin{2};
        Seq = VDJdata(:,SeqLoc);
        Header = VDJdata(:,HeaderLoc);
    end
    
    DataStruct = cell2struct([Header Seq],{'Header' 'Sequence'},2);
        
    DotLoc = regexp(FileName,'\.');
    SaveName = [FileName(1:DotLoc(end)-1) '.fa'];
    fastawrite([FilePath SaveName],DataStruct)  
    
    clear DataStruct VDJdata NewHeader SeqNum Header Seq
end


