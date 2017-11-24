%getGeneDatabase will get the gene database from a directory containing
%IMGT fasta files, labeled as IGXX_SPECIES.fa. 
%
%  DB = getGeneDatabase(Species)
%
%  DB = getGeneDatabase(FolderName)
%
%  Names = getGeneDatabase('getlist')
%
%  INPUT
%    Species: the species name, which should be the same as the folder name
%      store where this m file is. mfilepath/Species/
%    FolderName: Folder storing the IMGT fasta files (saved as
%      IGXX_SPECIES.fa). These fasta files should have the IMGT gap ("...")
%      in the V gene sequences to ensure all CDR and FWR regions are
%      aligned. Use this if using IMGT's fasta files, or if users want to
%      add sequences but with the gap information.
%
%    CsvFile: comma delimited file containing the IMGT genes without gap,
%      but marking where the gaps would have been and locations of CDR and
%      C and F/W anchor. This file is generated everytime someone runs this
%      using the FolderPath option instead.
%
%  OUTPUT
%    DB: Structure of database file containing the following fields:
%      FilePath     File path of this database
%      FileName     File name of this database csv file
%      MapHeader    1xN cell of column names
%      Vmap         V heavy gene DB
%      Dmap         D heavy gene DB
%      Jmap         J heavy gene DB
%      Vkmap        V kappa light gene DB
%      Jkmap        J kappa light gene DB
%      Vlmap        V lambda light gene DB
%      Jlmap        J lambda light gene DB
%
%  NOTE on map
%    Each map contains the following information:
%      'Seq': gene sequence without gap
%      'GeneName': gene name
%      'Function': functionality f, orf, p, [f], [p]
%      'Strain': mouse strain this gene was isolated from
%      'EntryNum': the Nth entry number in the original fasta file
%      'CDR1start': nt position where CDR1 starts
%      'CDR1end': nt position where CDR1 ends
%      'CDR2start': nt position where CDR2 starts
%      'CDR2end': nt position where CDR2 ends
%      'AnchorDist': number of nt from either 3' (or 5') side V (or J) to
%        the FIRST codon nt of the conserved 104C (or the 118F/W)
%      'GapInfo: location of the IMGT gaps that were removed, written as a
%        string N1-M1;N2-M2;...; where N is the nt position of the
%        ungapped sequence, and M is the number of gaps that follow AFTER
%        this position.
%
%  NOTE on CsvFile
%    This function will always check to see if there is a csv file in the
%    folder that is being summoned. If there is no csv file, it will
%    generate one based on the fasta files within that folder, and use it.
%
%  NOTE on gene locations
%    This function will look in the directory ./Database_Manager/[Species]
%    to load the sequences files stored in the csv file.

function DB = getGeneDatabase(varargin)
%Locate the database path
RootPath = findRoot();
DBPath = [RootPath 'Databases' filesep];
if ~exist(DBPath, 'dir')
    error('%s: Could not locate the "Databases" folder at [ %s ]', mfilename, DBPath);
end

%Determine the potential database folders
DBFolders = dir(DBPath);
DelLoc = zeros(length(DBFolders), 1, 'logical');
for j = 1:length(DelLoc)
    if ~DBFolders(j).isdir
        DelLoc(j) = 1;
    elseif strcmpi(DBFolders(j).name, '.') || strcmpi(DBFolders(j).name, '..')
            DelLoc(j) = 1;
    else %Check if there is a fa or csv file inside
        FastaFiles = dir([DBPath DBFolders(j).name filesep '*.fa']);
        CsvFiles = dir([DBPath DBFolders(j).name filesep '*.csv']);
        if isempty(FastaFiles) && isempty(CsvFiles)
            DelLoc(j) = 1;
        end
    end
end
DBFolders(DelLoc) = [];
DatabaseNames = struct2cell(DBFolders);
DatabaseNames = DatabaseNames(1, :)';

if ~isempty(varargin) && ischar(varargin{1}) && strcmpi(varargin{1}, 'getlist')
    DB = DatabaseNames;
    return;
end

if isempty(varargin) || isempty(varargin{1})
    SpeciesList = getGeneDatabase('getList');
    fprintf('What species is it?\n');
    dispList(SpeciesList);
    Attempt = 0;
    while 1
        Selection = input('Select option: ', 's');
        try 
            Selection = round(eval(Selection));
            if Selection > 0 && Selection <= length(SpeciesList)
                Species = SpeciesList{Selection};
                break;
            end
        catch
        end
        Attempt = Attempt + 1;
        if Attempt >= 5 
            error('%s: Did not choose correct option.', mfilename);
        end
    end
elseif ~isempty(varargin)
    Species = varargin{1};
else    
    error('%s: Input is incorrect.', mfilename);
end

%Determine which folder the user is specifying
DatabaseLoc = zeros(length(DBFolders), 1, 'logical');
for j = 1:length(DBFolders)
    if ~isempty(regexpi(DBFolders(j).name, Species, 'once'))
        DatabaseLoc(j) = 1;
    end
end
DatabaseIdx = find(DatabaseLoc);

%Make sure only 1 database can be selected
if isempty(DatabaseIdx) 
    error('%s: No database for this species found in [ %s ] .\n', mfilename, DBPath);
elseif length(DatabaseIdx) > 1
    dispList(DatabaseNames(DatabaseIdx));
    error('%s: Multiple databases found for this species. Be more specific, or use exact folder name\n.', mfilename);
end

DatabaseFolder = [DBPath DatabaseNames{DatabaseIdx} filesep];

%Get DB and display reference acknowledgement
DB = processIMGTfasta(DatabaseFolder);

if isempty(varargin) || ~any(contains(varargin, 'suppress', 'ignorecase', true))
    fprintf('  Germline gene databases were downloaded from http://www.imgt.org.\n');
    fprintf('  IMGT founder and director: Marie-Paule Lefranc, Montpellier, France\n');
    fprintf('\n');
end
