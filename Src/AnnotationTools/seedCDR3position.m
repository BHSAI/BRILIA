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

function VDJdata = seedCDR3position(VDJdata, Map, DB, X, Nleft, Nright, CheckSeqDir)
if contains(X, {'Vk', 'Vl'}, 'ignorecase', true)
    IsJ = 0;
    Chain = 'l';
    CIdx = 3; 
    SpecialSeed = {'TGT'};   %Conserved C
elseif contains(X, {'Jl', 'Jk'}, 'ignorecase', true)
    IsJ = 1;
    Chain = 'l';
    CIdx = 4;
    SpecialSeed = {'TTT'; 'TTC'};  %Conserved F for light chain
elseif strcmpi(X, 'V')
    IsJ = 0;
    Chain = 'h';
    CIdx = 3;
    SpecialSeed = {'TGT'};   %Conserved C
else
    IsJ = 1;
    Chain = 'h';
    CIdx = 4;   
    SpecialSeed = {'TGG'};   %Conserved W
end
if ~contains(Map.Chain, Chain, 'ignorecase', true); return; end

SeedPat = [sprintf('%s|', SpecialSeed{1:end-1}) SpecialSeed{end}];
SeqIdx = Map.([Chain 'Seq']);
CDR3Idx = Map.([Chain 'CDR3'])(CIdx);

Xseed = getGeneSeed(DB, X, Nleft, Nright, 'nt');
if isempty(Xseed); return; end

CheckRev = strcmpi(CheckSeqDir, 'y');

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

%Find the CDR3start or CDR3end locations. 
parfor j = 1:size(VDJdata, 1)
    Tdata = VDJdata(j, :);  
    if length(Tdata{SeqIdx}) <= (Nleft + Nright); continue; end
    [Score, StartAt] = alignSeqMEX(Tdata{SeqIdx}, Xseed, MissRate, Alphabet, ExactMatch, TrimSide, PenaltySide, PreferSide); 
    UnqPos = unique(StartAt(2, :) + Nleft + (StartAt(2, :) < 0));
    
    if CheckRev
        SeqR = seqrcomplement(Tdata{SeqIdx});
        [ScoreR, StartAtR] = alignSeqMEX(SeqR, Xseed, MissRate, Alphabet, ExactMatch, TrimSide, PenaltySide, PreferSide); 
        UnqPosR = unique(StartAtR(2, :) + Nleft + (StartAtR(2, :) < 0));
        
        if max(ScoreR(2, :)) > max(Score(2, :)) && length(UnqPosR) <= length(UnqPos) %Accept complement sequence if it's a higher score AND less positions
            Tdata{SeqIdx} = SeqR;
            UnqPos = UnqPosR;
        end
    end
    
    SpecPos = strfind(Tdata{SeqIdx}, SeedPat); %Include special CDR3 anchor locations
    CDR3Pos = unique([UnqPos(:); SpecPos(:)])';
    if IsJ; CDR3Pos = CDR3Pos +2; end %Need to include 2 nt of codon for J (since it's the end)
    ValidLoc = ~(CDR3Pos > length(Tdata{SeqIdx}))  |  (CDR3Pos < 1);
    CDR3Pos = CDR3Pos(ValidLoc);

    if ~isempty(CDR3Pos)
        Tdata{CDR3Idx} = CDR3Pos;
        VDJdata(j, :) = Tdata; %Add Tdata to sliced variable VDJdata
    end
end
