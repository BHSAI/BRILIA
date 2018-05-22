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
%    CellData: cell array (but must not have an embeded cell in a cell);
%    OutputFile: full name of the file to be saved to
%    Delimiter [',' ';' '\t']: Default is ','
%    'append': appends to the outfile instead of overwriting.
%
%  See also readDlmFile

function writeDlmFile(CellData, OutputFile, varargin)
if ismatrix(CellData) && ~iscell(CellData)
    CellData = num2cell(CellData);
end

%See if 'append' and Delimiter was specified
WriteType = 'w';
Delimiter = ',';
for j = 1:length(varargin)
    if ischar(varargin{j})
        if ~isempty(regexpi(varargin{j}, 'append', 'once'))
            WriteType = 'a';
        elseif ismember(varargin{j}, {',' ';' '\t'})
            Delimiter = varargin{j};
        else
            error('%s: Unrecognized input [ %s ].', mfilename, varargin{j});
        end
    end
end

%Convert all cell inputs as strings, and ensure no cell has the delimiter.
if strcmpi(Delimiter, '\t')
    ReplaceStr = ' ';
else
    ReplaceStr = Delimiter;
end
parfor j = 1:length(CellData(:))
    if isnumeric(CellData{j})
        if isnan(CellData{j}) 
            CellData{j} = 0;
        elseif isempty(CellData{j})
            CellData{j} = '';
        else
            CellData{j} = mat2str(CellData{j});
        end
    else
        CellData{j} = strrep(CellData{j}, ReplaceStr, '_');
    end
end

%Write the CellData as a formated file
TxtFormat = [repmat(['%s' Delimiter], 1, size(CellData, 2) - 1) '%s\n'];
[FID, MSG] = fopen(OutputFile, WriteType);
if FID < 0
    error('%s: Could not write to file [ %s ].\n  %s', mfilename, OutputFile, MSG);
end
try
    for j = 1:size(CellData, 1)
        fprintf(FID, TxtFormat, CellData{j, :});
    end
    fclose(FID);
catch
    fclose(FID); %Prevents memory leak
end
