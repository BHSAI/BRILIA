%readFasta will read a fasta file and return 2 cells containing the
%sequence name and the sequence itself. [Not used as MATLAB has the
%fastaread function]
%
%  [Header, Sequence] = readFasta(FileName);
%
%  [Header, Sequence] = readFasta();
%
%  INPUT
%    FileName: full file name of the fasta file. If empty, will ask user to
%    select the file.
%
%  OUTPUT
%    Header: Nx1 cell of sequence header information
%    Sequence: Nx1 cell of sequences

function [Header, Sequence] = readFasta(varargin)
SeqRange = [1 Inf];
for j = 1:length(varargin)-1
    if ischar(varargin{j}) && strcmpi(varargin{j},'seqrange') || strcmpi(varargin{j},'blockread')
        SeqRange = varargin{j+1};
        varargin(j:j+1) = [];
    end
end

%Get the input file name
if isempty(varargin)
    [FileName, FilePath] = uigetfile('*.txt;*.fa*','Open FASTA file');
    InputFileName = [FilePath FileName];
elseif ischar(varargin{1})
    InputFileName = varargin{1};
else
    error('readFasta: input file not specified correctly');
end

fprintf('%s: Fasta FileName: %s\n', mfilename, InputFileName);
fprintf('%s: SeqRange: %d to %d\n', mfilename, SeqRange(1), SeqRange(2));

%Determine number of sequences first
FID = fopen(InputFileName,'r');
if FID < 0
    error('%s: Could not open FASTA file, %s', mfilename, InputFileName);
end
NumLines = 0;
while 1
    TextLine= fgetl(FID);
    if ischar(TextLine) 
        if ~isempty(TextLine) && strcmpi(TextLine(1), '>')
            NumLines = NumLines+1;
        end
    else
        break
    end
end
fseek(FID,0,'bof');
fprintf('%s: Num of seq: %d\n', mfilename, NumLines);

if SeqRange(1) < 1; SeqRange(1) = 1; end
if SeqRange(2) > NumLines; SeqRange(2) = NumLines; end

%Start the matrix generation
Header = cell(diff(SeqRange) + 1,1);
Sequence = cell(diff(SeqRange) + 1,1);
j = 1;
CurLine = 0;
while CurLine < SeqRange(1)
    TextLine= fgetl(FID);
    if ischar(TextLine) 
        if ~isempty(TextLine) && strcmpi(TextLine(1), '>')
            CurLine = CurLine + 1;
        end
    else
        break
    end
end
if CurLine == 0
    return;
end
while CurLine <= SeqRange(2)
    if ~isempty(TextLine) && strcmp(TextLine(1),'>') %Get the header
        Header{j} = TextLine(2:end);
        SeqStr = '';
        TextLine = fgetl(FID);
        while TextLine(1) ~= '>' %Get the sequence
            SeqStr = cat(2,SeqStr,TextLine);
            TextLine = fgetl(FID);
            if isnumeric(TextLine); break; end;
            if isempty(TextLine); break; end;
        end
        Sequence{j} = strrep(SeqStr,' ','');
        j = j+1;
        CurLine = CurLine + 1;
    elseif isnumeric(TextLine)
        break;
    else        
        TextLine = fgetl(FID);
    end
end
fclose(FID);
