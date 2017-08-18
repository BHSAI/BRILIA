%writeGeneDatabaseToCsv will write a gene database to a csv file. The csv
%file is used mainly for debugging purposes to see if BRILIA is processing
%IMGT fasta files correctly and identifying the correct locations of CDR
%and FWR regions.
%
%  writeGeneDatabaseToCsv(DB, OutputFile)
%
%  writeGeneDatabaseToCsv(Species, OutputFile)
%
%  INPUT
%    DB: structured database output from processIMGTfasta
%    Species: folder name or species name of the database to print
%    OutputFile: full file name to write the csv file to
%  
%  OUTPUT
%    A csv file storing all no-gapped gene sequences used by BRILIA
%    processing.
%

function writeGeneDatabaseToCsv(Input, OutputFile)
if isstruct(Input)
    DB = Input;
elseif ischar(Input)
    DB = getGeneDatabase(Input);
end

%Determine maxmium number of entries (M) and data columns (N)
FieldNames = fieldnames(DB);
MapHeaderLoc = findCell(FieldNames, 'MapHeader');
MapLoc = findCell(FieldNames, 'map', 'MatchWord', 'Partial');
MapHeader = DB.(FieldNames{MapHeaderLoc});
MapNames = FieldNames(MapLoc);
    
%Determine data cell size for preallocation
M = 1; %Number of entries
for j = 1:length(MapNames)
    M = M + size(DB.(MapNames{j}), 1); 
end
N = size(MapHeader, 2); %Number of data columns

%Fill in the full table
AllData = cell(M, N);
AllData(1, :) = MapHeader;
S1 = 2;
for j = 1:length(MapNames)
    S2 = S1 + size(DB.(MapNames{j}), 1) - 1;
    AllData(S1:S2, :) = DB.(MapNames{j});
    S1 = S2 + 1;
end

%Ensure output file has a csv at the end
[FilePath, FileName, FileExt] = parseFileName(OutputFile, 'ignorefilecheck');
if ~strcmpi(FileExt, '.csv')
    FileName = [FileName '.csv'];
    FileName = strrep(FileName, '..', '.');
end
writeDlmFile(AllData, [FilePath FileName], ',');
