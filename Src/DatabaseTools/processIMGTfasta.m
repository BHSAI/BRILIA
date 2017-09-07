%processIMGTfasta will take the fasta germline VDJ gene sequences from 
%http://www.imgt.org/vquest/refseqh.html and convert them into usable
%database file for BRILIA. Below is the direction to use this.
%
%  1) Save reference sequences from http://www.IMGT.org/vquest/refseqh.html
%     into a fasta text file, where the file name starts with either:
%        IGHV_SpeciesName.fa
%        IGHD_SpeciesName.fa
%        IGHJ_SpeciesName.fa
%        IGKV_SpeciesName.fa
%        IGKJ_SpeciesName.fa
%        IGLV_SpeciesName.fa
%        IGLJ_SpeciesName.fa
%
%  2) Run processIMGTfasta, and select the folder where all those files are
%     stored.
%
%  OUTPUT
%     A csv files with the name [SpeciesName].csv storing all the processed
%       sequence file. This is what's loaded by getGeneDatabase.
%  
%  NOTE
%    There are 8 variables as follows:
%      Header: contains the header information of the ref seq maps
%      Vmap  = IGHV genes
%      Dmap  = IGHD genes
%      Jmap  = IGHJ genes
%      Vkmap = IGKV genes, kappa
%      Jkmap = IGKJ genes, kappa
%      Vlmap = IGLV genes, lambda
%      Jlmap = IGLJ genes, lambda
%
%    The Header contains the following information fields:
%      Seq: nucleotide sequence with IMGT's ... notation removed
%      GeneName: IMGT or custom gene name
%      Function: functionality of gene, P = pseudo, F = functional, ORF =
%        open reading frame.
%      Host: Host strain information, if available
%      CDR1start: Locations of CDR1 first nt
%      CDR1end: Locations of CDR1 last nt
%      CDR2start: Locations of CDR2 first nt
%      CDR2end: Locations of CDR2 last nt
%      AnchorDist: Distance from V gene 3' end, or J gene 5' end, to
%        the 1st nt of the codon of the constant residue 104C or 118 W/F.
%
%  NOTE
%    The the last column will ALWAYS be the anchor distance, due to other
%    function calling on the last column to be the anchor. Col 1 is ALWAYS
%    the sequence, and Col2 is ALWAYS the ref gene name.
%
%  DEVELOPER'S NOTE
%    The first and last column of a map matrix MUST be the Seq and
%    AnchorDist. This is because alignment algorithms assume this, to
%    prevent having to spend time searching for this.

function DB = processIMGTfasta(varargin)
%Find the folder where the database is stored
if isempty(varargin) || isempty(varargin{1})
    DatabaseFolder = uigetdir(cd, 'Select directory storing fasta files');
else
    DatabaseFolder = varargin{1};
end
if strcmpi(DatabaseFolder(end), filesep)
    DatabaseFolder = [DatabaseFolder filesep];
end

%Determine the species name based on the folder name
SlashLoc = regexpi(DatabaseFolder, filesep);
SpeciesName = DatabaseFolder(SlashLoc(end-1)+1:SlashLoc(end)-1);
SpeciesName = strrep(SpeciesName, ' ', '_');

%Find the GeneTable folder and file if it exists.
if exist([DatabaseFolder 'GeneTable'], 'dir') > 0
    GeneTableFiles = dir([DatabaseFolder 'GeneTable' filesep '*.csv']);
else
    GeneTableFiles = [];
end

%Order this
FileOrder = {'IGHV', 'IGHD', 'IGHJ', 'IGKV', 'IGKJ', 'IGLV', 'IGLJ'};
FastaFiles = dir([DatabaseFolder '*.fa*']);
DBlist = cell(7, 2); %fasta file, gene table file
for j = 1:length(FastaFiles)
    IGname = regexpi(FastaFiles(j).name, 'IG[HKL][VDJ]', 'match');
    if ~isempty(IGname)
        CellLoc = findCell(FileOrder, IGname, 'MatchCase', 'Any');
        if CellLoc(1) > 0
            DBlist{CellLoc, 1} = FastaFiles(j).name;
        end
    end
    
    if ~isempty(GeneTableFiles)
        IGname = regexpi(GeneTableFiles(j).name, 'IG[HKL][VDJ]', 'match');
        if ~isempty(IGname)
            IGname = upper(IGname{1});
            CellLoc = findCell(FileOrder, IGname, 'MatchCase', 'Any');
            if CellLoc(1) > 0
                DBlist{CellLoc, 2} = [DatabaseFolder 'GeneTable' filesep GeneTableFiles(j).name];
            end
        end
    end
end

%Determine the current fasta file check sum
CurCheckSum = 0;
for f = 1:size(DBlist, 1)
    if isempty(DBlist{f}) || ~ischar(DBlist{f, 1})
        continue; 
    end
    TempDir = dir(DBlist{f, 1});
    if ~isempty(TempDir)
        CurCheckSum = CurCheckSum + TempDir.bytes;
    end
end

%Set DB to empty if changes to database are detected
CsvFileName = [DatabaseFolder SpeciesName '.csv'];
MatFileName = [DatabaseFolder SpeciesName '.mat'];
DB = [];
FastaCheckSum = 0;
if exist(MatFileName, 'file')
    load(MatFileName);
    if FastaCheckSum ~= CurCheckSum
        DB = []; %Delete to alert this to redo processing
    else
        if ~exist(CsvFileName, 'file') && ~isempty(DB)
            writeGeneDatabaseToCsv(DB, CsvFileName);
        end
    end
end

%Process everything from scrap
if isempty(DB)
    fprintf('%s: Processing IMGT fasta files due to detected changes\n', mfilename);
    %--------------------------------------------------------------------------
    %V gene CDR start and end locations
    %Ref: http://www.imgt.org/IMGTScientificChart/Numbering/CDR1-IMGTgaps.html
    %Ref: http://www.imgt.org/IMGTScientificChart/Numbering/CDR2-IMGTgaps.html

    %CDR1 loc is the 27 to 38 AA residues
    CDR1start = 26*3 + 1;
    CDR1end = 38*3;

    %CDR2 loc is the 56 to 65 AA residues
    CDR2start = 55*3 + 1;
    CDR2end = 65*3;

    %CDR3 loc is the 105 AA residue, BUT, BRILIA uses the 104C anchor instead.
    CDR3anchor = 103*3 + 1;

    %--------------------------------------------------------------------------
    %Fill in the database information

    %Setup the DB
    MapHeader = {'Seq', 'GeneName', 'Function', 'Strain', 'EntryNum', 'GapInfo', 'CDR1start', 'CDR1end', 'CDR2start', 'CDR2end', 'AnchorDist'}; %Last 4 only applies for V genes
    DB.FilePath = DatabaseFolder;
    DB.MapHeader = MapHeader;
    MapNames = {'Vmap', 'Dmap', 'Jmap', 'Vkmap', 'Jkmap', 'Vlmap', 'Jlmap'};
    for j = 1:length(MapNames)
        DB.(MapNames{j}) = [];
    end

    %Setup the NT, AA, Name map.
    SeqLoc = findCell(MapHeader, 'Seq');
    GeneNameLoc = findCell(MapHeader, 'GeneName');
    EntryNumLoc = findCell(MapHeader, 'EntryNum');
    CDR3Loc = findCell(MapHeader, 'AnchorDist');
    for f = 1:size(DBlist, 1)
        if ~isempty(DBlist{f, 1})
            [SeqName, Seq] = fastaread(DBlist{f, 1});
            Xmap = cell(length(SeqName), length(MapHeader));

            %Load the gene table if it exists
            if ~isempty(DBlist{f, 2})
                GeneTable = readDlmFile(DBlist{f, 2}, 'Delimiter', ';');
            else
                GeneTable = [];
            end

            for j = 1:length(SeqName)
                %Extract seq info
                CurSeq = upper(Seq{j});               
                SeqInfo = regexpi(SeqName{j}, '\|', 'split');
                GeneName = SeqInfo{2};
                Function = SeqInfo{4};

                %Need to do this for the V's
                if max(f == [1, 4, 6]) > 0
                    %Get the default locations first.
                    CDR1s = CDR1start;
                    CDR1e = CDR1end;
                    CDR2s = CDR2start;
                    CDR2e = CDR2end;
                    CDR3a = CDR3anchor;

                    %FR1 can sometimes have 26a, 26b, etc. Need to shift CDRs
                    AddOn = 0;
                    while CDR1s <= length(CurSeq) && CurSeq(CDR1s) == '.'
                        AddOn = AddOn + 1;
                        CDR1s = CDR1s + 1;
                    end
                    if AddOn > 0 
                        CDR1e = CDR1e + AddOn;
                        CDR2s = CDR2s + AddOn;
                        CDR2e = CDR2e + AddOn;
                        CDR3a = CDR3a + AddOn;
                    end

                    %FR2 can sometimes have 55a, 55b, etc. Need to shift CDRs.
                    AddOn = 0;
                    while CDR2s <= length(CurSeq) && CurSeq(CDR2s) == '.'
                        AddOn = AddOn + 1;
                        CDR2s = CDR2s + 1;
                    end
                    if AddOn > 0
                        CDR2e = CDR2e + AddOn;
                        CDR3a = CDR3a + AddOn;
                    end

                    %FR3 can sometimes have gap. Need to shift CDRs.
                    AddOn = 0;
                    while CDR3a <= length(CurSeq) && CurSeq(CDR3a) == '.'
                        AddOn = AddOn + 1;
                        CDR3a = CDR3a + 1;
                    end

                    %Create a CDR tracking mat
                    TempSeq = zeros(size(Seq{j}));
                    TempSeq(CDR1s:CDR1e) = 1;
                    TempSeq(CDR2s:CDR2e) = 2;
                    TempSeq(CDR3a:end) = 3;

                    %Look for single ".", which are apparently wildcards...
                    OneDotLoc = regexp(CurSeq, '[^\.][\.]{1, 1}[^\.]') + 1;
                    TwoDotLoc = regexp(CurSeq, '[^\.][\.]{2, 2}[^\.]') + 1;
                    TwoDotLoc = cat(1, TwoDotLoc, TwoDotLoc+1);
                    TwoDotLoc = TwoDotLoc(:)';
                    CurSeq([OneDotLoc TwoDotLoc]) = 'X';
                    if ~isempty(OneDotLoc) || ~isempty(TwoDotLoc)
                        fprintf('%s has a < 3-gap segment\n', GeneName);
                    end

                    %Delete out the dots
                    [NoGapSeq, GapInfo] = removeIMGTgaps(CurSeq);
                    SeqDotLoc = CurSeq == '.';
                    TempSeq(SeqDotLoc) = [];
                    CurSeq = NoGapSeq;

                    %Determine new CDR locations
                    CDR1idx = find(TempSeq == 1);
                    CDR2idx = find(TempSeq == 2);
                    CDR3idx = find(TempSeq == 3);

                    %Determine CDR start and end locations
                    if ~isempty(CDR1idx)
                        CDR1s = CDR1idx(1);
                        CDR1e = CDR1idx(end);
                    end
                    if ~isempty(CDR2idx)
                        CDR2s = CDR2idx(1);
                        CDR2e = CDR2idx(end);
                    end
                    if ~isempty(CDR3idx)
                        CDR3anc = CDR3idx(1);
                    end

                    %Recalc CDR3anchor as 1st codon nt distance from 3' end
                    AnchorDist = length(CurSeq) - CDR3anc + 1;
                    if AnchorDist < 0
                        AnchorDist = 0;
                    end
                else
                    CDR1s = 0;
                    CDR1e = 0;
                    CDR2s = 0;
                    CDR2e = 0;
                    AnchorDist = 0;
                    GapInfo = '';
                end

                %Fill in the Strain data
                AllStrain = 'All';
                if ~isempty(GeneTable)
                    AllStrain = findStrainMatch(GeneName, GeneTable);
                end

                %Fill in the basic info in Xmap
                Xmap(j, :) = {CurSeq, GeneName, Function, AllStrain, j, GapInfo, CDR1s, CDR1e, CDR2s, CDR2e, AnchorDist};
            end

            %Need to adjust Xmap for J AnchorDist
            if max(f == [3, 5, 7]) > 0
                [Jalign, JconsSeq] = alignMultSeq(Seq');
                TempSeq = zeros(1, length(JconsSeq));

                %Set the search pattern depending on 118W or 118F
                if f == 3 %For IGHJ
                    Pattern = 'TGG'; %W
                else %For IGKJ
                    Pattern = 'TT[TC]GG[ACGT]'; %FG
                end

                %Find where the conserved residue appears most
                for j = 1:length(SeqName)
                    MatchLoc = regexpi(Jalign{j}, Pattern);
                    TempSeq(MatchLoc) = TempSeq(MatchLoc) + 1;
                end
                AnchorLoc = find(TempSeq == max(TempSeq));
                AnchorLoc = AnchorLoc(1);

                %Establish the CDR3 anchor point
                CDR3anchors = zeros(length(Seq), 1);
                for j = 1:length(SeqName)
                    CurSeq = Jalign{j};
                    CurSeq(AnchorLoc) = '!'; %Just a marker;
                    CurSeq = strrep(CurSeq, '-', '');
                    CDR3anchors(j) = find(CurSeq == '!');
                end

                %Fill in the Xmap anchor
                Xmap(:, CDR3Loc) = num2cell(CDR3anchors);
            end

            %Need to adjust Xmap for D fwd and rev
            if f == 2
                TempMap = cell(size(Xmap, 1)*2, size(Xmap, 2));
                TempMap(1:2:end, :) = Xmap;
                for j = 1:size(Xmap, 1)
                    Xmap{j, SeqLoc} = seqrcomplement(Xmap{j, 1});
                    Xmap{j, GeneNameLoc} = ['r' Xmap{j, 2}];
                end
                TempMap(2:2:end, :) = Xmap;
                Xmap = TempMap;
                Xmap(:, EntryNumLoc) = num2cell(1:size(Xmap,1));
                clear TempMap;
            end
            DB.(MapNames{f}) = Xmap;
        end
    end

    %Create the csv format for easy viewing
    writeGeneDatabaseToCsv(DB, CsvFileName);

    %Save this DB as a .mat file for faster loading
    FastaCheckSum = CurCheckSum;
    save(MatFileName, 'DB', 'FastaCheckSum');
end
