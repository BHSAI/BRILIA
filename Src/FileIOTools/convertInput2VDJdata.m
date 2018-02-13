%convertInput2VDJdata will assess what file input is being dealt with and
%then extract the relevant information required to make the VDJdata format
%file used by BRILIA. This initializes the VDJdata cell. 
%
%  [VDJdata, VDJheader] = convertInput2VDJdata()
%
%  [VDJdata, VDJheader] = convertInput2VDJdata(FullFileName)
%
%  [VDJdata, VDJheader] = convertInput2VDJdata(FullFileName, 'FileType', FileType, 'Delimiter', Delimiter)
%
%  INPUT
%    FullFileName: Full name of input file. If empty, will ask users. 
%    FileType ['fasta', 'fastq', 'delimited']: Specifying
%      FileType prevents erroneous determination of file type.
%    Delimiter [';' ',' '\t']: Needed only for delimited file type
%
%  OUTPUT
%    VDJdata: main BRILIA data cell
%    VDJheader: main BRILIA header cell (edit at for Headers_BRILIA.csv)

function [VDJdata, VDJheader, varargout] = convertInput2VDJdata(varargin)
%Parse the inputs
P = inputParser;
addOptional(P, 'FullFileName', '', @(x) ischar(x) || isempty(x));
addParameter(P, 'Chain', 'H', @(x) ismember({upper(x)}, {'H', 'L', 'HL'}));
addParameter(P, 'FileType', '', @(x) ismember({lower(x)}, {'', 'fasta', 'fastq', 'delimited'}));
addParameter(P, 'Delimiter', '', @(x) ismember({lower(x)}, {';', ',', '\t', ''})); 
addParameter(P, 'SeqRange', [1 Inf], @(x) isnumeric(x));
parse(P, varargin{:});
P = P.Results;
P.SeqRange = double(P.SeqRange); %For some reason, this is uint64. So change to double.

%Determine the file type here
if isempty(P.FullFileName)
    [InFileName, InFilePath] = uigetfile('*.fa*;*.*sv', 'Select the input sequence file', 'MultiSelect', 'off');
    if isempty(InFileName)
        error('%s: No file was selected.', mfilename);
    end
    P.FullFileName = [InFilePath InFileName];
end
[InFilePath, InFileName, InFileExt] = parseFileName(P.FullFileName);    

%Check to see if there is a file type override
if isempty(P.FileType)
    if strcmpi(InFileExt, '.fa') || strcmpi(InFileExt, '.fasta')
        P.FileType = 'fasta';
    elseif strcmpi(InFileExt, '.fastq')
        P.FileType = 'fastq';
    else
        P.FileType = 'delimited';
    end
end

%==========================================================================
%Begin assembling the VDJdata structure

%Determine what data the input files have
if strcmpi(P.FileType, 'fasta')
    %Open the fasta file and convert to cell
    try
        [SeqName, SeqData] = readFasta(P.FullFileName, 'SeqRange', P.SeqRange);
    catch
        error('%s: Could not read fasta file ''%s''.', mfilename, P.FullFileName);
    end
    if ischar(SeqName)
        SeqName = {SeqName};
        SeqData = {SeqData};
    end
    InputData = [SeqName(:) SeqData(:)];
    InSeqNameLoc = 1;
    InSeqLoc = 2;
    InTemplateLoc = 0;
elseif strcmpi(P.FileType, 'fastq')
    %Open the fasta file and convert to cell
    [SeqName, SeqData] = fastqread(P.FullFileName, 'blockread', P.SeqRange);
    if ischar(SeqName)
        SeqName = {SeqName};
        SeqData = {SeqData};
    end
    InputData = [SeqName(:) SeqData(:)];        
    InSeqNameLoc = 1;
    InSeqLoc = 2;
    InTemplateLoc = 0;
elseif strcmpi(P.FileType, 'delimited')
    InputHeader = readDlmFile(P.FullFileName, 'delimiter', P.Delimiter, 'LineRange', [1 1]); %For getting the header only
    InputData = readDlmFile(P.FullFileName, 'delimiter', P.Delimiter, 'LineRange', P.SeqRange + 1); %For LineRange, we assume 1st line is header, hence get the next one.

    %Note: Only delimited files can take in paired sequences
    InTemplateLoc = find(contains(InputHeader, {'Template', 'TempCount', 'TempCt', 'Copy'}, 'IgnoreCase', true));
    InSeqNameLoc = find(contains(InputHeader, 'SeqName', 'IgnoreCase', true));
    InHSeqLoc = find(ismember(lower(strrep(InputHeader, '-', '')), 'hseq'));
    InLSeqLoc = find(ismember(lower(strrep(InputHeader, '-', '')), 'lseq'));
    InXSeqLoc = find(ismember(lower(strrep(InputHeader, '-', '')),  'seq'));
    
    if isempty(InHSeqLoc) && isempty(InLSeqLoc) && isempty(InXSeqLoc)
        error('%s: Could not find the "Seq" or "H-Seq" or "L-Seq" column header in the file "%s".', mfilename, P.FullFileName);
    end
    
    if isempty(InHSeqLoc) && isempty(InLSeqLoc)
        InSeqLoc = InXSeqLoc;
    elseif isempty(InHSeqLoc) && ~isempty(InLSeqLoc)
        InSeqLoc = InLSeqLoc;
        P.Chain = 'L';
    elseif ~isempty(InHSeqLoc) && isempty(InLSeqLoc) 
        InSeqLoc = InHSeqLoc;
        P.Chain = 'H';
    else
        InSeqLoc = [InHSeqLoc InLSeqLoc];
        P.Chain = 'HL';
    end
    
    if length(InSeqLoc) > 3
        error('%s: Too many "Seq" or "H-Seq" or "L-Seq" in the column header in the file "%s".', mfilename, P.FullFileName);
    end
end
        
%Create the VDJdata default matrix
[VDJdata, VDJheader] = getBlankDataTable(size(InputData, 1), P.Chain);
[H, L, ~] = getAllHeaderVar(VDJheader);

%Fill in VDJdata
VDJdata(:, H.TemplateLoc) = num2cell(ones(size(InputData, 1), 1)); %Always initialize TempCt column with 1.
VDJdata(:, H.SeqNumLoc) = num2cell([1:size(InputData, 1)] + P.SeqRange(1) - 1); %Always assign a unique numbering order for VDJdata
VDJdata(:, H.GrpNumLoc) = num2cell([1:size(InputData, 1)] + P.SeqRange(1) - 1); %Always assign a unique numbering order for VDJdata
VDJdata(:, H.SeqNameLoc) = InputData(:, InSeqNameLoc);
if InTemplateLoc > 0 %Update template loc if possible.
    for j = 1:size(VDJdata, 1)
        if isnumeric(InputData{j, InTemplateLoc})
            VDJdata(j, H.TemplateLoc) = InputData(j, InTemplateLoc);
        elseif ischar(InputData{j, InTemplateLoc}) %Check if it is a number
            if min(isstrprop(InputData{j, InTemplateLoc}, 'digit')) == 1
                VDJdata{j, H.TemplateLoc} = convStr2Num(InputData{j, InTemplateLoc});
            end
        end
    end
end

%Copy the sequences
if strcmpi(P.Chain, 'H')
    VDJdata(:, H.SeqLoc) = InputData(:, InSeqLoc);
elseif strcmpi(P.Chain, 'L')
    VDJdata(:, L.SeqLoc) = InputData(:, InSeqLoc);
else
    VDJdata(:, H.SeqLoc) = InputData(:, InSeqLoc(1));
    VDJdata(:, L.SeqLoc) = InputData(:, InSeqLoc(2));
end

if nargout >= 3
    varargout{1} = InFileName;
    if nargout >= 4
        varargout{2} = InFilePath;
    end
end
