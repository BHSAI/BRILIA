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

function VDJdata = seedCDR3position(VDJdata, VDJheader, DB, X, Nleft, Nright, CheckSeqDir)
CheckSeqDir = upper(CheckSeqDir(1)); %Ensure proper format

%Determine the chain and segment
Xparse = regexpi(X, ',\s*|;', 'split');
if strcmpi(Xparse{1}, 'V') || strcmpi(Xparse{1}, 'J')
    if length(Xparse) > 1
        error('%s: Cannot have multiple heavy chains. Check X variable.', mfilename);
    end
    Chain = 'H';
else
    Chain = 'L';
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
[H, L, ~] = getAllHeaderVar(VDJheader);
if Chain == 'H' %Heavy chain
    SeqLoc = H.SeqLoc;
    if SeqLoc == 0; return; end %Invalid 
    if Segment == 'V'
        SpecialSeed = {'TGT'};      %conserved C
        CDR3colLoc = H.CDR3Loc(3);  %Where in VDJdata to store the location
    else
        SpecialSeed = {'TGG'};      %Conserved W
        CDR3colLoc = H.CDR3Loc(4);  %Where in VDJdata to store the location
    end
else %Light chain
    SeqLoc = L.SeqLoc;
    if SeqLoc == 0; return; end %Invalid 
    if Segment == 'V'
        SpecialSeed = {'TGT'};      %conserved C
        CDR3colLoc = L.CDR3Loc(3);  %Where in VDJdata to store the location
    else
        SpecialSeed = {'TTT'; 'TTC'};   %Conserved F for light chain
        CDR3colLoc =  L.CDR3Loc(4);     %Where in VDJdata to store the location
    end
end
RepPat = repmat('%s|', 1, length(SpecialSeed));
RepPat(end) = [];
SeedPat = sprintf(RepPat, SpecialSeed{:});  %Special seed for V and J, which is the 'TGT' and the 'TGG' or 'TT[TC]'

%Setup the alignSeq input structure for faster alignment
P.MissRate = 0;
P.Alphabet = 'nt';
P.CheckSeq = 'yes';
P.DiagIdx = [];
P.ExactMatch = 'no';
P.PreferSide = 'none';
if Segment == 'V'
    P.TrimSide = 'right';   %will trim poorly matched regions in this side
    P.PreferSide = 'left';  %if tied alignments, will favor one towards this side of seed and seq
    P.PenaltySide = 'left'; %unused nts of seed in this side will be unfavorable
elseif Segment == 'J'
    P.TrimSide = 'left';    %will trim poorly matched regions in this side
    P.PreferSide = 'right'; %if tied alignments, will favor one towards this side of seed and seq
    P.PenaltySide = 'right';%unused nts of seed in this side will be unfavorable
end

%Find the CDR3start or CDR3end locations. Flip seq too if needed.
parfor j = 1:size(VDJdata, 1)
    %Extract a "slice" of VDJdata for parfor
    Tdata = VDJdata(j, :);  
    Seq = Tdata{SeqLoc};
    if length(Seq) < Nleft + Nright + 1; %If sequence is too short, skip
        continue; 
    end 
    
    %Do seed alignments, forward sense
    AlignScores = zeros(size(Xseed, 1), 3);
    for x = 1:size(Xseed, 1)
       [Score, ~, StartAt, MatchAt] = alignSeq(Xseed{x}, Seq, P); %Subtle hint: alignSeq 1st input should be shorter than 2nd input for speed
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
            [Score, ~, StartAt, MatchAt] = alignSeq(Xseed{x}, SeqNTR, P);
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
