%seedCDR3position will find "potential" locations for 104C first or 118W/F
%last nt codon within in a sequence. This function could also check both
%forward and reverse direction, flipping the sequence direction with the
%highest seed alignment score and lowest variability in seed location
%(since incorrect direction will lead to many bad alignments and
%positions).
%
%  VDJdata = seedCDR3position(VDJdata,VDJheader,DB,X,Nleft,Nright,CheckSeqDir)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    VDJheader: main BRILIA header cell
%    DB: Gene database structure (getCurrentDatabase.m)
%    X ['V' 'Vk' 'Vl' 'J' 'Jk' 'Jl']: Specificy what database gene is used.
%      Can use combinations such as 'Vl,Vk'.
%    Nleft: num of nts left of CDR3 codon 1st nt Anchor for align seed. 
%    Nright: num of nts right of CDR3 codon 1st nt Anchor for align seed.
%    CheckSeqDir ['n' or 'y']: no or yes for checking both sequence
%      directions. If 'y', Will take complements of sequences if needed.
%
%  OUTPUT
%    VDJdata: processed BRILIA data cell
%
%  NOTE
%    This will fill VDJdata with the CDR3 start and end positions with
%    single or multiple integer values. The findVDJmatch will use these
%    seeded locations to speed up the alignment process, whereas
%    updateVDJdata will select the final values.
%
%    Cannot seed a heavy and light chain V or J together. Ex, X = 'V, Vk'
%    will not work.
%
%  See also findVDJmatch, findVJmatch, findGeneMatch

function VDJdata = seedCDR3position(VDJdata, Map, DB, X, Nleft, Nright, CheckSeqDir)
CheckSeqDir = upper(CheckSeqDir(1)); %Ensure proper format

if contains(Map.Chain, 'H', 'IgnoreCase', true) && any(contains(X, {'k', 'l'}, 'IgnoreCase', true))
    return;
elseif contains(Map.Chain, 'L', 'IgnoreCase', true) && ~any(contains(X, {'k', 'l'}, 'IgnoreCase', true))
    return;
end

%Determine the chain and segment
Xparse = regexpi(X, ',\s*|;', 'split');
if ( strcmpi(Xparse{1}, 'V') || strcmpi(Xparse{1}, 'J') ) && length(Xparse) > 1
    error('%s: Cannot have multiple heavy chains. Check X variable.', mfilename);
end
Segment = upper(X(1));

%Extract the seed sequences
Xseed = {''};
for j = 1:length(Xparse)
    Xseed = cat(1, Xseed, getGeneSeed(DB, Xparse{j}, Nleft, Nright, 'nt'));
end
Xseed = unique(Xseed);
if isempty(Xseed{1})
    Xseed(1) = [];
end
if isempty(Xseed) %missing all database files or valid seed sequences
    return;
end

%Get header locations since parfor can't handle broadcast variables
if strcmpi(Map.Chain, 'H') %Heavy chain
    SeqLoc = Map.hSeq;
    if SeqLoc == 0; return; end %Invalid 
    if Segment == 'V'
        SpecialSeed = {'TGT'};      %conserved C
        CDR3colLoc = Map.hCDR3(3);  %Where in VDJdata to store the location
    else
        SpecialSeed = {'TGG'};      %Conserved W
        CDR3colLoc = Map.hCDR3(4);  %Where in VDJdata to store the location
    end
elseif strcmpi(Map.Chain, 'L')  %Light chain
    SeqLoc = Map.lSeq;
    if SeqLoc == 0; return; end %Invalid 
    if Segment == 'V'
        SpecialSeed = {'TGT'};      %conserved C
        CDR3colLoc = Map.lCDR3(3);  %Where in VDJdata to store the location
    else
        SpecialSeed = {'TTT'; 'TTC'};  %Conserved F for light chain
        CDR3colLoc = Map.lCDR3(4);     %Where in VDJdata to store the location
    end
else %Can't determine H or L side.
    return; 
end
SeedPat = sprintf('%s|', SpecialSeed{:});  %Special seed for V and J, which is the 'TGT' and the 'TGG' or 'TT[TC]'
SeedPat(end) = [];

%Setup the input for alignSeqMEX.
MissRate   = 0;
Alphabet   = 'n'; %nt
ExactMatch = 'n'; 
if X == 'V'
    TrimSide    = 'r'; %right
    PenaltySide = 'l'; %left
    PreferSide  = 'l'; %left
elseif X == 'J'
    TrimSide    = 'l'; %left
    PenaltySide = 'r'; %right
    PreferSide  = 'r'; %right
elseif X == 'D'
    TrimSide    = 'b'; %both
    PenaltySide = 'n'; %none
    PreferSide  = 'n'; %none
else
    TrimSide    = 'n'; %none
    PenaltySide = 'n'; %none
    PreferSide  = 'n'; %none
end

%Find the CDR3start or CDR3end locations. Flip seq too if needed.
parfor j = 1:size(VDJdata, 1)
    Tdata = VDJdata(j, :);  
    Seq = Tdata{SeqLoc};
    if length(Seq) < (Nleft + Nright + 1); continue; end %Skip short seqs
    
    %Do seed alignments, forward sense
    AlignScores = zeros(size(Xseed, 1), 3);
    for x = 1:size(Xseed, 1)
       [Score, StartAt, MatchAt,a] = alignSeqMEX(Xseed{x}, Seq, MissRate, Alphabet, ExactMatch, TrimSide, PenaltySide, PreferSide); 
       if StartAt(2) > 0
            AnchorLoc = Nleft - StartAt(2) + 2; %Comes from (Nleft + 1) - (StartAt(2) - 1) = PositionAnchorSeed - MissingSeqCount
       else
            AnchorLoc = Nleft + abs(StartAt(2)) + 1;
       end       
       AlignScores(x, :) = [Score(1)/(diff(MatchAt)+1) Score(2) AnchorLoc];
    end

    %Do seed alignments on reverse complement sequence
    if CheckSeqDir == 'Y'
        AlignScoresR = zeros(size(Xseed, 1), 3);
        SeqNTR = seqrcomplement(Seq);
        for x = 1:size(Xseed, 1)
            [Score, StartAt, MatchAt] = alignSeqMEX(Xseed{x}, SeqNTR, MissRate, Alphabet, ExactMatch, TrimSide, PenaltySide, PreferSide); 
            if StartAt(2) > 0
                AnchorLoc = Nleft - StartAt(2) + 2; %Comes from (Nleft + 1) - (StartAt(2) - 1) = PositionAnchorSeed - MissingSeqCount
            else
                AnchorLoc = Nleft + abs(StartAt(2)) + 1;
            end       
            AlignScoresR(x, :) = [Score(1)/(diff(MatchAt)+1) Score(2) AnchorLoc];
        end

        %Accept complement sequence if it's a higher score AND less positions
        if max(AlignScoresR(:, 2)) > max(AlignScores(:, 2))
            FwdPos = unique(AlignScores(:, 3));  %Check how many unique anchor positions there are
            RevPos = unique(AlignScoresR(:, 3));
            if length(RevPos) <= length(FwdPos) %Means seeds are converging more on reverse direction          
                Tdata{1, SeqLoc} = SeqNTR;
                AlignScores = AlignScoresR;
            end
        end
    end
    
    %Finalize all potential CDR3 anchor locations, including the special
    %conserved codon
    SpecPos = regexp(Seq, SeedPat);
    CDR3Pos = unique([AlignScores(:, 3); SpecPos(:)]);
    if Segment == 'J' %Need to include 2 nt of codon for J (since it's the end)
        CDR3Pos = CDR3Pos + 2; 
    end
    InvalidPos = (CDR3Pos > length(Seq))  |  (CDR3Pos < 1);
    CDR3Pos(InvalidPos) = [];

    %Update to VDJdata ONLY if there is a seed
    if ~isempty(CDR3Pos)
        Tdata{1, CDR3colLoc} = unique([Tdata{1, CDR3colLoc}(:); CDR3Pos(:)])';
        VDJdata(j, :) = Tdata; %Add Tdata to sliced variable VDJdata
    end
end
