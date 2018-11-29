%getGeneDatabase will get the gene database from a directory containing
%IMGT fasta files, labeled as IGXX_SPECIES.fa. 
%
%  DB = getGeneDatabase(SpeciesFolderName)
%
%  Names = getGeneDatabase('getlist')
%
%  INPUT
%    SpeciesFolderName: the name of the species, which is the same as the
%    Databases folder that is created where the BRILIA.exe file is. These
%    Database\SpeciesFolderName folders contains IMGT fasta files that are
%    saved as IGHV.fa, IGHD.fa, IGHJ.fa, IGKV.fa. etc. These fasta files
%    should have the IMGT gap ("...") in the V gene sequences to ensure all 
%    CDR and FWR regions are aligned. 
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
%  NOTE on using custom gene databases
%    This function will look in [Exe_Root]/Databases/[Species] for the
%    database folders. These folders should be named by their species name.
%    To use a custom database, create a new folder within the Databases
%    folder, and then add the fasta files. The files MUST start with IGHV,
%    IGHJ, IGHD, IGKV, IGKL, IGLV, IGLJ so that BRILIA knows how to process
%    these fasta files. There alsmot must be a complete VDJ set for heavy
%    chain, and VJ set for light chain.
%
%    Example: 
%    1) BRILIA.exe is located at:
%        C:\User\BRILIA
%    2) Custom mouse genes is added as: 
%        C:\User\BRILIA\Databases\mouse1
%    3) Heavy chain fasta files with IMGT gap notations are added as:
%        C:\User\BRILIA\Databases\mouse1\IGHV.fa
%                                       \IGHD.fa
%                                       \IGHJ.fa
%    4) Run BRILIA as: 
%        C:\User\BRILIA.exe my_seq_file.fastq Species mouse1 Chain H
%
%  NOTE on map
%    Each DB.Xmap contains the following cell array columns in order:
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
function [DB, FiltOption] = getGeneDatabase(varargin)
%Locate the databases
RootPath = findExeDir(); %Will find either EXE or code root dir
DBPath = fullfile(RootPath, 'Databases');
if isdeployed
    TempDBPath = fullfile(findRoot(), 'Databases'); %findRoot will find the temporary code root dir
    if ~exist(DBPath, 'dir')
        [Success, Msg] = copyfile(TempDBPath, DBPath);
        assert(Success, '%s: Could not copy file "%s" to "%s".\n %s', mfilename, TempDBPath, DBPath, Msg);
    else %Need to figure out what is NOT there
        [DirToAdd, DirNameToAdd] = dir2(TempDBPath, 'dir');
        [~, DirNameExist] = dir2(DBPath, 'dir');
        Idx = find(~ismember(DirNameToAdd, DirNameExist));
        for q = 1:length(Idx)
            copyfile(DirToAdd{Idx(q)}, fullfile(DBPath, DirNameToAdd{Idx(q)}));
        end
    end
end
assert(exist(DBPath, 'dir') > 0, '%s: Could not locate the "Databases" folder at "%s".', mfilename, DBPath);

[DBFolders, SpeciesList] = dir2(DBPath, 'dir');

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