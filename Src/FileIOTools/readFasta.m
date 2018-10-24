%readFasta will read a fasta file. Similar to the MATLAB fastaread.m, but
%can handle files with empty lines and odd header names.
%
%  S = readFasta();
%
%  [Header, Seq] = readFasta();
%
%  [...] = readFasta(FileName);
%
%  [...] = readFasta(FileName, 'SeqRange', Range)
%
%  INPUT
%    FileName: file name of the fasta file. If empty, will ask user to
%      select the file.
%    Range: 2-element integer vector specifying first and last sequence
%      to read.
%
%  OUTPUT
%    Header: Nx1 cell of sequence header information
%    Seq: Nx1 cell of sequences
%    S: Nx1 structure containing fields 'Header' and 'Seq'

function varargout = readFasta(varargin)
%Determine the ith to jth sequence range to read
SeqRange = [1 Inf]; %Set 1 to Inf as default. The Inf will be changed later.
OptionIdx = find(strcmpi(varargin, 'seqrange') | strcmpi(varargin, 'blockread'));
if ~isempty(OptionIdx)
    SeqRange = varargin{OptionIdx(1)+1};
    if numel(SeqRange) == 1
        SeqRange = repelem(SeqRange, 1, 2);
    end
    SeqRange(SeqRange < 1) = 1; 
    varargin(OptionIdx(1):OptionIdx(1)+1) = [];
end

%Check if user specified a file name
if isempty(varargin)
    [FileName, FilePath] = uigetfile('*.txt;*.fa*', 'Open FASTA file');
    InputFileName = fullfile(FilePath, FileName);
elseif exist(varargin{1}, 'file') || exist(fullfile(pwd, varargin{1}), 'file')
    InputFileName = varargin{1};
else
    error('%s: input file not specified correctly.', mfilename);
end

%Read the file (Note: This will read all sequence text. This could be slow
%if you have a large sequence file but only want a few sequences.
[FID, MSG] = fopen(InputFileName, 'r');
assert(FID > 0, '%s: Could not open FASTA file "%s".\n  %s', mfilename, InputFileName, MSG);
TXT = textscan(FID, '%s', 'delimiter', '\n');
TXT = TXT{1};
fclose(FID);

%Identify where the ">" are, which are the headers
Idx = [find(contains(TXT, '>')); length(TXT)+1];  %header index. Added "length(TXT)" to simplify for loop later
SeqRange(2) = min([SeqRange(2) length(Idx)-1]); %fix SeqRange(2)
SeqCount = diff(SeqRange) + 1;

if nargout == 1 %Output to a structure
    S(SeqCount) = struct('Header', '', 'Sequence', '');
    k = 1;
    for j = SeqRange(1):SeqRange(2)
        S(k).Header = strrep(TXT{Idx(j)}(find(TXT{Idx(j)} == '>', 1)+1:end), '"', '');
        S(k).Sequence = strrep(strcat(TXT{Idx(j)+1:Idx(j+1)-1}), ' ', '');
        k = k+1;
    end
    varargout{1} = S;
elseif nargout == 2 %Output to cells
    Header = cell(SeqCount, 1);
    Sequence = cell(SeqCount, 1);
    k = 1;
    for j = SeqRange(1):SeqRange(2)
        Header{k} = strrep(TXT{Idx(j)}(find(TXT{Idx(j)} == '>', 1)+1:end), '"', '');
        Sequence{k} = strrep(strcat(TXT{Idx(j)+1:Idx(j+1)-1}), ' ', '');
        k = k+1;
    end
    varargout{1} = Header;
    varargout{2} = Sequence;
else
    error('%s: Too many outputs. Maximum is 2.', mfilename);
end