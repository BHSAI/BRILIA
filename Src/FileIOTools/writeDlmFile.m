%writeDlmFile will take a cell array and convert it into a string array, 
%and then save it into a delimited file. If the delimiter is the same as a
%character in the string array, will automatically replace the string
%character with a "|" character to ensure delimiter number is not ruined.
%
%  writeDlmFile(CellData, OutputFile, Delimiter)
%
%  writeDlmFile(CellData, OutputFile, Delimiter, 'append')
%
%  INPUT
%    CellData: cell array (but must not have a cell in a cell);
%    OutputFile: full name of the file to be saved to
%    Delimiter [',' ';' '\t']: Default is ','
%    'append': appends to the outfile instead of overwriting.
%
%  See also readDlmFile

function writeDlmFile(CellData, OutputFile, varargin)
if ~iscell(CellData)
    CellData = num2cell(CellData);
end

WriteType = 'w';
Delimiter = ',';
if ~isempty(varargin)   
    AppendLoc = strcmpi(varargin, 'append');
    if any(AppendLoc)
        WriteType = 'a';
        varargin = varargin(~AppendLoc);
    end

    DelimLoc = strcmpi(varargin, {',', ';', '\t'});
    if any(DelimLoc)
        Delimiter = varargin{find(DelimLoc, 1)};
    end
end

CellData(cellfun('isempty', CellData)) = {''};
NumLoc = cellfun('isclass', CellData, 'double');
StrLoc = ~NumLoc;
CellData(NumLoc) = cellfun(@mat2str, CellData(NumLoc), 'un', 0); 
if contains(Delimiter, '\t') %To \t from becoming tab
    CellData(StrLoc) = strrep(CellData(StrLoc), '\t', '\\t');
else %Prevent string with delimiter from messing up delimiter setup
    CellData(StrLoc) = strrep(CellData(StrLoc), Delimiter, '|');
end

[FID, MSG] = fopen(OutputFile, WriteType);
assert(FID > 0, '%s: Could not create/write to file "%s".\n  %s', mfilename, OutputFile, MSG);
Fmt = [repmat(['%s' Delimiter], 1, size(CellData, 2)-1) '%s\n'];
for j = 1:size(CellData, 1)
    fprintf(FID, Fmt, CellData{j, :});
end
fclose(FID);