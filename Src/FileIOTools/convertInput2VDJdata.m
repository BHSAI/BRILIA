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
%     MinQuality  * '2'                    Phred Score (ASCII Base = 33) for P_error = 0.01995. Only for fastq files.
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
addParameter(P, 'Chain',     'H',     @(x) ismember({upper(x)}, {'H', 'L', 'HL', 'LH'}));
addParameter(P, 'FileType',  '',      @(x) ismember({lower(x)}, {'', 'fasta', 'fastq', 'delimited'}));
addParameter(P, 'Delimiter', '',      @(x) ismember({lower(x)}, {'', ';', ',', '\t'})); 
addParameter(P, 'SeqRange',  [1 Inf], @(x) isnumeric(x) && isa(x, 'double'));
addParameter(P, 'MinQuality', '2',    @(x) ischar(x) || isnumeric(x)); %ASCII_BASE=33, '2' = P_error 0.01995
parse(P, varargin{:});
P = P.Results;

FileName   = P.FileName;
Chain      = regexprep(upper(P.Chain), 'LH', 'HL');
FileType   = P.FileType;
Delimiter  = P.Delimiter;
SeqRange   = P.SeqRange;
MinQuality = P.MinQuality;
if isnumeric(MinQuality) %make sure it's a character for Phred Score ASCII base 33
    MinQuality = sprintf('%d', MinQuality);
end

%Get the file name
if isempty(FileName)
    [InFileName, InFilePath] = uigetfile('*.fa*;*.*sv', 'Select the input sequence file', 'MultiSelect', 'off');
    assert(ischar(InFileName), '%s: No file was selected.', mfilename);
    FileName = fullfile(InFilePath, InFileName);
end
[InFilePath, InFileName, InFileExt] = parseFileName(FileName);    

%Determine the FileType
if isempty(FileType)
    if strcmpi(InFileExt, '.fa') || strcmpi(InFileExt, '.fasta')
        FileType = 'fasta';
    elseif strcmpi(InFileExt, '.fastq')
        FileType = 'fastq';
    else
        FileType = 'delimited';
    end
end

%==========================================================================
%Begin assembling the VDJdata structure

%Determine what data the input files have
Template = {1};
switch lower(FileType)
    case 'fasta'
        [SeqName, SeqData] = readFasta(FileName, 'SeqRange', SeqRange);

    case 'fastq'
        [SeqName, SeqData, Quality] = fastqread(FileName, 'blockread', SeqRange);

        if ischar(SeqData) %In case the fastq file has only 1 sequence, output would be char.
            SeqName = {SeqName};
            SeqData = {SeqData};
            Quality = {Quality};            
        end
        
        %Convert low-quality nt reads as 'N'
        MinQuality = uint8(MinQuality);
        parfor j = 1:length(SeqName)
            SeqData{j}(uint8(Quality{j}) < MinQuality) = 'N'; 
        end

    case 'delimited'
        InHeader = readDlmFile(FileName, 'delimiter', Delimiter, 'LineRange', 1); %For getting the header only
        InHeader = formatStrSame(InHeader); %For string matching below
        InData   = readDlmFile(FileName, 'delimiter', Delimiter, 'LineRange', SeqRange + 1); %For LineRange, we assume 1st line is header, hence get the next one.

        %Note: Only delimited files can take in paired sequences
        InTemplateIdx = find(contains(InHeader, {'template', 'tempcount', 'tempct', 'copy'}), 1);
        InSeqNameIdx  = find(contains(InHeader, 'seqname'), 1);
        InHSeqIdx = find(ismember(InHeader, 'hseq'), 1);
        InLSeqIdx = find(ismember(InHeader, 'lseq'), 1);
        InXSeqIdx = find(ismember(InHeader,  'seq'), 1);

        if isempty(InHSeqIdx) && isempty(InLSeqIdx) && isempty(InXSeqIdx)
            error('%s: Could not find the "Seq" or "H-Seq" or "L-Seq" column header in the file "%s".', mfilename, FileName);
        end

        if isempty(InHSeqIdx) && isempty(InLSeqIdx) %Assume whatever Seq is the Chain's Seq
            InSeqIdx = InXSeqIdx;
        elseif isempty(InLSeqIdx) && ~isempty(InHSeqIdx) %Only H Chain exist. Use this.
            InSeqIdx = InHSeqIdx;
            if Chain ~= 'H'
                warning('%s: User set Chain = ''L'', but only ''H-Seq'' exists. Setting Chain to ''H''.', mfilename);
                Chain = 'H';
            end
        elseif isempty(InHSeqIdx) && ~isempty(InLSeqIdx) %Only L Chain exist. Use this.
            InSeqIdx = InLSeqIdx;
            if Chain ~= 'L'
                warning('%s: User set Chain = ''H'', but only ''L-Seq'' exists. Setting Chain to ''L''.', mfilename);
                Chain = 'L';
            end
        else
            InSeqIdx = [InHSeqIdx InLSeqIdx]; %Both chains exist. Choose the user-specified chain.
            if numel(Chain) < 2
                InSeqIdx = InSeqIdx('HL' == Chain);
            end
        end
        
        if ~isempty(InTemplateIdx)
            Template = cellfun(@convStr2NumMEX, InData(:, InTemplateIdx), 'un', 0);
        end
        SeqName = InData(:, InSeqNameIdx);
        SeqData = InData(:, InSeqIdx);
end
        
%Create the VDJdata default matrix
[VDJdata, VDJheader] = getBlankDataTable(size(SeqData, 1), Chain);
Map = getVDJmapper(VDJheader);

VDJdata(:, Map.SeqName) = SeqName;
VDJdata(:, Map.Template) = Template;
VDJdata(:, Map.SeqNum) = num2cell(1:size(SeqData, 1));
VDJdata(:, Map.GrpNum) = num2cell(1:size(SeqData, 1));
for c = 1:numel(Map.Chain)
    C = lower(Map.Chain(c));
    VDJdata(:, Map.([C 'Seq'])) = SeqData(:, c);    
end

%Format string to be same for matching purposes. 
%EDIT_NOTE: Edit this code if the string matching criteria changes
function Str = formatStrSame(Str)
Str = lower(regexprep(Str, '[^a-zA-Z0-9\|]', ''));