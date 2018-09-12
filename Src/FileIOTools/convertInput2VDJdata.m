%convertInput2VDJdata will assess what file input is being dealt with and
%then extract the relevant information required to make the VDJdata format
%file used by BRILIA. This initializes the VDJdata cell. 
%
%  [VDJdata, VDJheader, FileName, FilePath, Map] = convertInput2VDJdata()
%
%  [VDJdata, VDJheader] = convertInput2VDJdata(FileName)
%
%  [VDJdata, VDJheader] = convertInput2VDJdata(FileName, 'FileType', FileType, 'Delimiter', Delimiter)
%
%  INPUT
%    FileName: Full name of input file. If empty, will ask users. 
%
%     Param       Value (* = default)      Details
%     ----------- ------------------------ --------------------------------
%     FileType    * ''                     Autodetect sequence input file type
%                   'fasta'                Assume fasta files
%                   'fastq'                Assume fastq files
%                   'delimited'            Assume delmited file
%     Chain       * H                      Heavy chain
%                   L                      Light chain
%                   HL                     Heavy and Light chains
%     Delimiter   * ''                     Autodetect delimiter. Needed only for delimited file type.
%                   ','                    CSV
%                   ';'                    SSV
%                   '\t'                   TSV
%     SeqRange    * [1,Inf]                Process all sequences 
%                   #                      Process only the #th sequence
%                   [M,N]                  Process Mth to Nth seqeunce (include brackets "[]" , "," , and NO SPACE)
%     MinQuality  * '2'                    Pred Score (ASCII Base = 33) for P_error = 0.01995. Only for fastq files.
%                                          %https://www.drive5.com/usearch/manual/quality_score.html for info.
%
%  OUTPUT
%    VDJdata: main BRILIA data cell
%    VDJheader: main BRILIA header cell (obtained from DataHeaderInfo.csv)
%    FileName: file name without hte path
%    FilePath: file path
%    Map: structure of BRILIA data index 

function [VDJdata, VDJheader, InFileName, InFilePath, Map] = convertInput2VDJdata(varargin)
P = inputParser;
addOptional( P, 'FileName',  '',      @(x) isempty(x) || (ischar(x) && exist(x, 'file')));
addParameter(P, 'Chain',     'h',     @(x) ismember({lower(x)}, {'h', 'l', 'hl', 'lh'}));
addParameter(P, 'FileType',  '',      @(x) ismember({lower(x)}, {'', 'fasta', 'fastq', 'delimited'}));
addParameter(P, 'Delimiter', '',      @(x) ismember({lower(x)}, {'', ';', ',', '\t'})); 
addParameter(P, 'SeqRange',  [1 Inf], @(x) isnumeric(x) && isa(x, 'double'));
addParameter(P, 'MinQuality', '2',    @(x) ischar(x) || isnumeric(x)); %ASCII_BASE=33, '2' = P_error 0.01995
parse(P, varargin{:});
P = P.Results;
P.Chain = strrep(upper(P.Chain), 'LH', 'HL');

if isnumeric(P.MinQuality)
    P.MinQuality = sprintf('%d', P.MinQuality);
end

if isempty(P.FileName)
    [InFileName, InFilePath] = uigetfile('*.fa*;*.*sv', 'Select the input sequence file', 'MultiSelect', 'off');
    assert(ischar(InFileName), '%s: No file was selected.', mfilename);
    P.FileName = fullfile(InFilePath, InFileName);
end

[InFilePath, InFileName, InFileExt] = parseFileName(P.FileName);    

%Determine the File Type
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
        [SeqName, SeqData] = readFasta(P.FileName, 'SeqRange', P.SeqRange);
    catch
        error('%s: Could not read fasta file "%s".', mfilename, P.FileName);
    end
    InputData = [SeqName(:) SeqData(:)];
    InSeqNameLoc = 1;
    InSeqLoc = 2;
    InTemplateLoc = 0;
elseif strcmpi(P.FileType, 'fastq')
    %Open the fasta file and convert to cell
    if P.MinQuality == 0
        [SeqName, SeqData] = fastqread(P.FileName, 'blockread', P.SeqRange);
        if ischar(SeqName)
            SeqName = {SeqName};
            SeqData = {SeqData};
        end
    else
        [SeqName, SeqData, Quality] = fastqread(P.FileName, 'blockread', P.SeqRange);
        if ischar(Quality)
            SeqName = {SeqName};
            SeqData = {SeqData};
            Quality = {Quality};            
        end
        MinQuality = uint8(P.MinQuality);
        parfor j = 1:length(SeqName)
            BadLoc = uint8(Quality{j}) < MinQuality;
            SeqData{j}(BadLoc) = 'N';
        end
    end
    InputData = [SeqName(:) SeqData(:)];        
    InSeqNameLoc = 1;
    InSeqLoc = 2;
    InTemplateLoc = 0;
elseif strcmpi(P.FileType, 'delimited')
    InputHeader = readDlmFile(P.FileName, 'delimiter', P.Delimiter, 'LineRange', [1 1]); %For getting the header only
    InputData = readDlmFile(P.FileName, 'delimiter', P.Delimiter, 'LineRange', P.SeqRange + 1); %For LineRange, we assume 1st line is header, hence get the next one.

    %Note: Only delimited files can take in paired sequences
    InTemplateLoc = find(contains(InputHeader, {'Template', 'TempCount', 'TempCt', 'Copy'}, 'IgnoreCase', true));
    InSeqNameLoc = find(contains(InputHeader, 'SeqName', 'IgnoreCase', true));
    InHSeqLoc = find(ismember(lower(strrep(InputHeader, '-', '')), 'hseq'));
    InLSeqLoc = find(ismember(lower(strrep(InputHeader, '-', '')), 'lseq'));
    InXSeqLoc = find(ismember(lower(strrep(InputHeader, '-', '')),  'seq'));
    
    if isempty(InHSeqLoc) && isempty(InLSeqLoc) && isempty(InXSeqLoc)
        error('%s: Could not find the "Seq" or "H-Seq" or "L-Seq" column header in the file "%s".', mfilename, P.FileName);
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
        error('%s: Too many "Seq" or "H-Seq" or "L-Seq" in the column header in the file "%s".', mfilename, P.FileName);
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

Map = getVDJmapper(VDJheader);
