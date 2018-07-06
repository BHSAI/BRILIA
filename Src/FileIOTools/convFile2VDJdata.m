%convFile2VDJdata will assess what file input is being dealt with and
%then extract the relevant information required to make the VDJdata format
%file used by BRILIA.
%
%  VDJdata = convFile2VDJdata()
%
%  VDJdata = convFile2VDJdata(FullFileName)
%
%  VDJdata = convFile2VDJdata(FullFileName, 'FileType', FileType, 'Delimiter', Delimiter)
%
%  INPUT
%    FullFileName: Full name of input file. If empty, will ask users. 
%    FileType ['fasta', 'fastq', 'delimited']: Specifying
%      FileType prevents erroneous determination of file type.
%    Delimiter [';' ',' '\t']: Needed only for delimited file type
%
%  OUTPUT
%    VDJdata: BRILIA data structure

function VDJdata = convFile2VDJdata(FileName, varargin)
Chain     = parseInputVar('Chain', 'H', @(x)any(strcmpi(x, {'H', 'L', 'HL'})), varargin{:});
FileType  = parseInputVar('FileType', '', @(x)any(strcmpi(x, {'', 'fasta', 'fastq', 'delimited'})), varargin{:});
Delimiter = parseInputVar('Delimiter', '', @(x)any(strcmpi(x, {';', ',', '\t', ''})), varargin{:}); 
SeqRange  = parseInputVar('SeqRange', [1 Inf], @(x)isnumeric(x), varargin{:});

%Determine the file type here
if nargin == 0 || ~(exist(FileName, 'file') || exist(fullfile(pwd, FileName), 'file'))
    FileName = '';
end
if isempty(FileName)
    FullFileName = openFileDialog('*.fa*;*.*sv', 'Select the input sequence file', 'MultiSelect', 'off');
    FileName = FullFileName{1};
    if isempty(FileName)
        error('%s: No file was selected.', mfilename);
    end
end
[~, ~, FileExt] = parseFileName(FileName);

%Check to see if there is a file type override
if isempty(FileType)
    if any(strcmpi({'.fa', '.fasta'}, FileExt))
        FileType = 'fasta';
    elseif strcmpi(FileExt, '.fastq')
        FileType = 'fastq';
    else
        FileType = 'delimited';
    end
end

%Begin assembling the VDJdata structure
switch FileType
    case 'fasta'
        try
            [SeqName, SeqData] = readFasta(FileName, 'SeqRange', SeqRange);
        catch
            error('%s: Could not read fasta file "%s".', mfilename, FileName);
        end
    case 'fastq'
        try
            [SeqName, SeqData] = fastqread(FileName, 'blockread', SeqRange);
        catch
            error('%s: Could not read fastq file "%s".', mfilename, FileName);
        end
        if ischar(SeqName)
            SeqName = {SeqName};
            SeqData = {SeqData};
        end
    otherwise
        InputData = readDlmFile(FileName, 'delimiter', Delimiter, 'LineRange', SeqRange + 1); %For LineRange, we assume 1st line is header, hence get the next one.
        InputHeader = readDlmFile(FileName, 'delimiter', Delimiter, 'LineRange', [1 1]); %For getting the header only
        InputHeader = lower(InputHeader);
        
        HSeqLoc = ismember(lower(InputHeader), {'hseq', 'h-seq'});
        LSeqLoc = ismember(lower(InputHeader), {'lseq', 'l-seq'});
        ASeqLoc = ismember(lower(InputHeader), 'seq');
        CountLoc = endsWith(InputHeader, {'copynumber', 'templatecount'}, 'ignorecase', true);
        SeqNameLoc = endsWith(InputHeader, {'name', 'id'}, 'ignorecase', true);
        if ~any(HSeqLoc) && ~any(LSeqLoc)
            if any(ASeqLoc)
                HSeqLoc = ASeqLoc;
            else
                error('%s: Could not identify the sequence column. Label the column header as "hSeq" or "lSeq".', mfilename);
            end
        end
        
        
        
        %Note: Only delimited files can take in paired sequences
        Htemp = getHeavyHeaderVar(InputHeader);
        Ltemp = getLightHeaderVar(InputHeader);
        InSeqNameLoc = max([Htemp.SeqNameLoc Ltemp.SeqNameLoc]);       
        if InSeqNameLoc == 0; InSeqNameLoc = 1; end %Assume.

        InTemplateLoc = max([Htemp.TemplateLoc Ltemp.TemplateLoc]);
        if InTemplateLoc == 0 && size(InputData, 2) >= 4
            InTemplateLoc = 4; %Assume this
        end 

        if strcmpi(Chain, 'HL')        
            InSeqLoc = Htemp.SeqLoc;
            if InSeqLoc == 0; InSeqLoc = 2; end %Assume this

            InSeqLocL = Ltemp.SeqLoc;
            if InSeqLocL == 0; InSeqLoc = 3; end %Assume this
        else
            InSeqLoc = max([Htemp.SeqLoc Ltemp.SeqLoc]);
            if InSeqLoc == 0; InSeqLoc = 2; end %Assume this
        end        
end
        
%Create the VDJdata default matrix
[VDJdata, VDJheader] = getBlankDataTable(size(InputData, 1), Chain);
[H, L, ~] = getAllHeaderVar(VDJheader);

%Fill in VDJdata
VDJdata(:, H.TemplateLoc) = num2cell(ones(size(InputData, 1), 1)); %Always initialize TempCt column with 1.
VDJdata(:, H.SeqNumLoc) = num2cell([1:size(InputData, 1)] + SeqRange(1) - 1); %Always assign a unique numbering order for VDJdata
VDJdata(:, H.GrpNumLoc) = num2cell([1:size(InputData, 1)] + SeqRange(1) - 1); %Always assign a unique numbering order for VDJdata
VDJdata(:, H.SeqNameLoc) = InputData(:, InSeqNameLoc);
if InTemplateLoc > 0 %Update template loc if possible.
    for j = 1:size(VDJdata, 1)
        if isnumeric(InputData{j, InTemplateLoc})
            VDJdata(j, H.TemplateLoc) = InputData(j, InTemplateLoc);
        elseif ischar(InputData{j, InTemplateLoc}) %Check if it is a number
            if min(isstrprop(InputData{j, InTemplateLoc}, 'digit')) == 1
                VDJdata{j, H.TemplateLoc} = convStr2NumMEX(InputData{j, InTemplateLoc});
            end
        end
    end
end

%Copy the sequences
if strcmpi(Chain, 'H')
    VDJdata(:, H.SeqLoc) = InputData(:, InSeqLoc);
elseif strcmpi(Chain, 'L')
    VDJdata(:, L.SeqLoc) = InputData(:, InSeqLoc);
else
    VDJdata(:, H.SeqLoc) = InputData(:, InSeqLoc);
    VDJdata(:, L.SeqLoc) = InputData(:, InSeqLocL);
end