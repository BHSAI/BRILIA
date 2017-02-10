%convertInput2VDJdata will assess what file input is being dealt with and
%then extract the relevant information required to make the VDJdata format
%file used by BRILIA. This initializes the VDJdata cell. 
%
%  [VDJdata,VDJheader] = convertInput2VDJdata()
%
%  [VDJdata,VDJheader] = convertInput2VDJdata(FullFileName)
%
%  [VDJdata,VDJheader] = convertInput2VDJdata(FullFileName,'FileType',FileType,'Delimiter',Delimiter)
%
%  INPUT
%    FullFileName: Full name of input file. If empty, will ask users to
%      look for it. 
%    FileType ['fasta','fastq','excel','delimited']: Specifying
%      FileType prevents erroneous determination of file type.
%    Delimiter [';' , ',' , '\t']: Needed only for delimited file type
%
%  OUTPUT
%    VDJdata: main cell matrix used by BRILIA to store data
%    VDJheader: Name of each data column of VDJdata. To modify this, look
%      for Headers_BRILIA.csv

function [VDJdata,VDJheader,varargout] = convertInput2VDJdata(varargin)
P = inputParser;
addOptional(P,'FullFileName','',@(x) ischar(x) || isempty(x));
addParameter(P,'FileType','',@(x)ismember({lower(x)},{'','fasta','fastq','excel','delimited'}));
addParameter(P,'Delimiter',';',@(x)ismember({lower(x)},{';' ',' '\t'})); 
parse(P,varargin{:});
FullFileName = P.Results.FullFileName;
FileType = P.Results.FileType;
Delimiter = P.Results.Delimiter;

%Determine the file type here
if isempty(FullFileName)
    [InFileName, InFilePath] = uigetfile('*.fa*;*.xls*;*.csv','Select the input sequence file','MultiSelect','off');
    if isempty(InFileName);
        return
    end
    FullFileName = [InFilePath InFileName];
end
[InFilePath, InFileName, InFileExt] = parseFileName(FullFileName);    

%Check to see if there is a file type override
if isempty(FileType)
    if strcmpi(InFileExt,'.fa') || strcmpi(InFileExt,'.fasta')
        FileType = 'fasta';
    elseif strcmpi(InFileExt,'.fastq')
        FileType = 'fastq';
    elseif strcmpi(InFileExt(1:4),'.xls')
        FileType = 'excel';
    elseif strcmpi(InFileExt,'.csv') || strcmpi(InFileExt,'.tsv')
        FileType = 'delimited';
    else
        error('FileType cannot be determined. Make sure file ext is correct');
    end
end

%==========================================================================
%Begin assembling the VDJdata structure

%Determine what data the input files have
if strcmpi(FileType,'fasta')
    %Open the fasta file and convert to cell
    [SeqName,SeqData] = fastaread(FullFileName);
    InputData = [SeqName(:) SeqData(:)];
    InSeqNameLoc = 1;
    InSeqLoc = 2;
    InTemplateLoc = 0;
elseif strcmpi(FileType,'fastq')
    %Open the fasta file and convert to cell
    [SeqName,SeqData] = fastqread(FullFileName);
    InputData = [SeqName(:) SeqData(:)];        
    InSeqNameLoc = 1;
    InSeqLoc = 2;
    InTemplateLoc = 0;
elseif strcmpi(FileType,'excel') || strcmpi(FileType,'delimited')
    if strcmpi(FileType,'excel')
        [~,~,InputData] = xlsread(FullFileName);
    else
        InputData = readDlmFile(FullFileName,'delimiter',Delimiter);
    end
    InputHeader = InputData(1,:); %Assume 1st row is headers
    InputData(1,:) = [];

    %See if this the same as VDJdata header format
    InSeqNameLoc = findCell(InputHeader,'SeqName');
    InSeqLoc = findCell(InputHeader,{'nucleotide','Seq'});
    InTemplateLoc = findCell(InputHeader,{'count (templates)','TemplateCount'});
    if InSeqNameLoc == 0 || InSeqLoc == 0 || InTemplateLoc == 0 %Not VDJdata, assume standard format.
        InSeqNameLoc = 1;
        InSeqLoc = 2;
        if size(InputData,2) >= 3
            InTemplateLoc = 3;
        else
            InTemplateLoc = 0;
        end
    end        
end
        
%Create the VDJdata default matrix
HeaderData = readDlmFile('Headers_BRILIA.csv','Delimiter',';'); %Obtain the VDJdata header info for output format
VDJheader = HeaderData(2:end,1)';
H = getHeaderVar(VDJheader);
VDJdata = cell(size(InputData,1),length(VDJheader)); 
VDJdata(:,H.TemplateLoc) = num2cell(ones(size(InputData,1),1)); %Always initialize TempCt column with 1.
VDJdata(:,H.SeqNumLoc) = num2cell(1:size(InputData,1)); %Always assign a unique numbering order for VDJdata
VDJdata(:,H.GrpNumLoc) = num2cell(1:size(InputData,1)); %Always assign a unique numbering order for VDJdata
        
%Fill in the information that you have
VDJdata(:,H.SeqNameLoc) = InputData(:,InSeqNameLoc);
VDJdata(:,H.SeqLoc) = InputData(:,InSeqLoc);
if InTemplateLoc > 0 %Update template loc if possible.
    for j = 1:size(VDJdata,1)
        if isnumeric(InputData{j,InTemplateLoc})
            VDJdata(j,H.TemplateLoc) = InputData(j,InTemplateLoc);
        elseif ischar(InputData{j,InTemplateLoc}) %Check if it is a number
            if min(isstrprop(InputData{j,InTemplateLoc},'digit')) == 1
                VDJdata{j,H.TemplateLoc} = eval(InputData{j,InTemplateLoc});
            end
        end
    end
end

%Outputs
if nargout >= 3
    varargout{1} = InFileName;
    if nargout >= 4
        varargout{2} = InFilePath;
    end
end
