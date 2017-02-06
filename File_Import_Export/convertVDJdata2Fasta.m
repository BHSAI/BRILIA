%convertVDJdata2Fasta will ask  users to select multiple files of the
%VDJdata format, and then convert them to corresponding FASTA files for
%processing. By default, sequence will be named according to the "SeqNum"
%column if  it exists. if not, it will do a simple numerical ordering.

function convertVDJdata2Fasta(varargin)
[FileNames, FilePath] = uigetfile('*.xlsx;*.csv','Select the sequence data files','Multiselect','on');
if ischar(FileNames)
    FileNames = {FileNames};
end

for j = 1:length(FileNames)
    FileName = FileNames{j};
    [VDJdata,VDJheader] = openSeqData([FilePath FileName]);    
    H = getHeaderVar(VDJheader);
    
    %Assign seq header name if none is inputted
    if isempty(varargin)
        if H.SeqLoc == 0
            error('Error: Could not find the nucleotide column');
        else
            Seq = VDJdata(:,H.SeqLoc);
        end

        if H.SeqNumLoc == 0
            disp('Warning: Could not find the SeqNum column - making own numbering');
            SeqNum = num2cell([1:size(VDJdata,1)]');
        else
            SeqNum = VDJdata(:,H.SeqNumLoc);
        end

        Header = cell(size(VDJdata,1),1);
        for k = 1:size(Header,1)
            Header{k} = sprintf('Seq_%d',SeqNum{k});
        end
    else
        HeaderLoc = varargin{1};
        H.SeqLoc = varargin{2};
        Seq = VDJdata(:,H.SeqLoc);
        Header = VDJdata(:,HeaderLoc);
    end
    
    DataStruct = cell2struct([Header Seq],{'Header' 'Sequence'},2);
        
    DotLoc = regexp(FileName,'\.');
    SaveName = [FileName(1:DotLoc(end)-1) '.fa'];
    fastawrite([FilePath SaveName],DataStruct)  
    
    clear DataStruct VDJdata VDJheader SeqNum Header Seq
end
