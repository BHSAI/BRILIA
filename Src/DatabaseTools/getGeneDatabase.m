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
%    FiltOpt: Structure showing what filters were used for the database
%      Species      Species filter used
%      Strain       Strain filter used (only for mouse)
%      Vgene        V gene function filter used (ex: orf,f,p)
%      Dgene        D gene direction (fwd vs inv) filter used
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

function [DB, FiltOption] = getGeneDatabase(varargin)
%Locate the databases
RootPath = findRoot();
DBPath = fullfile(RootPath, 'Databases'); [RootPath 'Databases' filesep];
if ~exist(DBPath, 'dir')
    error('%s: Could not locate the "Databases" folder at "%s"', mfilename, DBPath);
end

DBFolders = dir(DBPath);
DBFolders = DBFolders(~(ismember({DBFolders.name}, {'.', '..'}) | ~[DBFolders.isdir]));
DBFolders = fullfile({DBFolders.folder}, {DBFolders.name});
KeepLoc = ~cellfun(@(y) ~any(arrayfun(@(x) endsWith(x.name, {'.fa', '.csv'}, 'ignorecase', true), dir(y))), DBFolders);
DBFolders = DBFolders(KeepLoc);
SpeciesList = cellfun(@(x) x(find(x == filesep, 1, 'last')+1:end), DBFolders, 'un', 0);

%Return for list search only
if ~isempty(varargin) && ischar(varargin{1}) && any(strcmpi(varargin{1}, {'getlist', 'list'}))
    DB = SpeciesList(:);
    return
end

%Select the species
if isempty(varargin) || isempty(varargin{1})
    SpeciesIdx = chooseFromList(SpeciesList, 'Attempt', 5, 'Default', [], 'Message', 'What species is it?');
else
    SpeciesLoc = strcmpi(SpeciesList, varargin{1});
    if ~any(SpeciesLoc) %Maybe try contains instead for nearest match
        SpeciesLoc = contains(SpeciesList, varargin{1}, 'ignorecase', true);
    end
    SpeciesIdx = find(SpeciesLoc);
end
if isempty(SpeciesIdx)
    error('%s: Did not choose correct Species option.', mfilename);
elseif length(SpeciesIdx) > 1
    error('%s: Species name is not specific enough. Multiple database possible for "%s".', mfilename, varargin{1});
end

%Process request to get database
DB = processIMGTfasta(DBFolders{SpeciesIdx});
if nargin > 1
    [DB, FiltOption] = filterGeneDatabase(DB, varargin{2:end});
end
FiltOption.Species = SpeciesList{SpeciesIdx};