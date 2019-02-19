%convertInput2VDJdata will assess what file input is being dealt with and
%then extract the relevant information required to make the VDJdata format
%file used by BRILIA. This initializes the VDJdata cell. For delimited
%files, the required header names are: 'Seq', 'HSeq', or 'LSeq'. Optional
%headers are 'SeqName' and 'Template'.
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
P.addOptional('FileName',  '',       @(x) isempty(x) || (ischar(x) && ~isempty(dir(x))));
P.addParameter('Chain',     'H',     @(x) any(strcmpi(x, {'H', 'L', 'HL', 'LH'})));
P.addParameter('FileType',  '',      @(x) any(strcmpi(x, {'', 'fasta', 'fastq', 'delimited'})));
P.addParameter('Delimiter', '',      @(x) any(strcmpi(x, {'', ';', ',', '\t'}))); 
P.addParameter('SeqRange',  [1 Inf], @(x) isnumeric(x) && numel(x) <= 2);
P.addParameter('MinQuality', '2',    @(x) ischar(x) || isnumeric(x)); %ASCII_BASE=33, '2' = P_error 0.01995
P.parse(varargin{:});

FileName   = P.Results.FileName;
Chain      = P.Results.Chain;
FileType   = P.Results.FileType;
Delimiter  = P.Results.Delimiter;
SeqRange   = P.Results.SeqRange;
MinQuality = P.Results.MinQuality;

%Make sure Chain is all caps and HL, not LH
Chain = regexprep(upper(Chain), 'LH', 'HL');
%Make sure MinQuality is character for Phred Score ASCII base 33
if isnumeric(MinQuality) 
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
    if any(strcmpi(InFileExt, {'.fa', 'fasta'}))
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
        SeqCount = numel(SeqName);
        Template = repelem({1}, SeqCount, 1);

    case 'fastq'
        [SeqName, SeqData, Quality] = fastqread(FileName, 'blockread', SeqRange);
        if ischar(SeqData) %In case the fastq file has only 1 sequence, output would be char.
            SeqName = {SeqName};
            SeqData = {SeqData};
            Quality = {Quality};
        end
        SeqCount = numel(SeqName);
        Template = repelem({1}, SeqCount, 1);
        
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
        
        %Determine the sequence nucleotide loctaion in delimited file
        InHSeqIdx = find(ismember(InHeader, 'hseq'), 1);
        InLSeqIdx = find(ismember(InHeader, 'lseq'), 1);
        InXSeqIdx = find(ismember(InHeader,  'seq'), 1);
        if isempty(InHSeqIdx) && isempty(InLSeqIdx) && isempty(InXSeqIdx)
            error('%s: Could not find the "Seq", "H-Seq", or "L-Seq" column header in "%s".', mfilename, FileName);
        end
        if isempty(InHSeqIdx) && isempty(InLSeqIdx) %Assume whatever Seq is the Chain's Seq
            InSeqIdx = InXSeqIdx;
        elseif isempty(InLSeqIdx) && ~isempty(InHSeqIdx) %Only H Chain exist. Use this.
            InSeqIdx = InHSeqIdx;
            if Chain ~= 'H' %Maybe user specified wrong Chain. Give warning but do not change, as it could be intentional.
                warning('%s: User set Chain = ''L'', but only ''H-Seq'' exists. Setting Chain to ''H''.', mfilename);
                Chain = 'H';
            end
        elseif isempty(InHSeqIdx) && ~isempty(InLSeqIdx) %Only L Chain exist. Use this.
            InSeqIdx = InLSeqIdx;
            if Chain ~= 'L' %Maybe user specified wrong Chain. Give warning but do not change, as it could be intentional. 
                warning('%s: User set Chain = ''H'', but only ''L-Seq'' exists. Setting Chain to ''L''.', mfilename);
                Chain = 'L';
            end
        else
            InSeqIdx = [InHSeqIdx InLSeqIdx]; %Both chains exist. Choose the user-specified chain.
            if numel(Chain) < 2
                InSeqIdx = InSeqIdx('HL' == Chain);
            end
        end
        
        %Determine the sequence name, if provided        
        InSeqNameIdx  = find(contains(InHeader, 'seqname'), 1);
        if ~isempty(InSeqNameIdx)
            SeqName = InData(:, InSeqNameIdx);
        else
            SeqName = cellfun(@(x) sprintf('%d', x), num2cell(1:size(InData, 1))', 'un', 0);
        end
        SeqData = InData(:, InSeqIdx);
        SeqCount = numel(SeqName);
        
        %Determine the sequence template counts, if provided
        InTemplateIdx = find(contains(InHeader, {'template', 'tempcount', 'tempct', 'copy'}), 1);
        if ~isempty(InTemplateIdx)
            Template = cellfun(@convStr2NumMEX, InData(:, InTemplateIdx), 'un', 0);
        else
            Template = repelem({1}, SeqCount, 1);
        end
end
        
%Create the VDJdata default matrix
[VDJdata, VDJheader] = getBlankDataTable(SeqCount, Chain);
Map = getVDJmapper(VDJheader);

VDJdata(:, Map.SeqName)  = SeqName;
VDJdata(:, Map.Template) = Template;
VDJdata(:, Map.SeqNum)   = num2cell(1:SeqCount);
VDJdata(:, Map.GrpNum)   = num2cell(1:SeqCount);
for c = 1:numel(Map.Chain)
    C = lower(Map.Chain(c));
    VDJdata(:, Map.([C 'Seq'])) = SeqData(:, c);    
end

%Format string to be same for matching purposes. 
%EDIT_NOTE: Edit this code if the string matching criteria changes
function Str = formatStrSame(Str)
Str = lower(regexprep(Str, '[^a-zA-Z0-9\|]', ''));