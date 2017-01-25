%seedCDR3position will find potential locations where 104C's 1st or 118W's
%last nt codon could reside in a sequence. This could also check both
%forward and reverse direction, flipping the sequence direction with the
%highest seed alignment score and lowest variability in seed location
%(since incorrect direction will lead to many bad alignments and
%positions).
%
%  VDJdata = seedCDR3position(VDJdata,NewHeader,Xmap,X,Nleft,Nright,CheckSeqDir)
%
%  INPUT
%    Xmap [Vmap or Dmap]: Database for V or J germline gene
%    X ['V' or 'J']: Specificy what database gene is used
%    Nleft: number of nts left of CDR3start's or CDR3end's, 1st nt of the
%      codon to use as alignment seed.
%    Nright: number of nts right of 104C to use as alignment seed.
%    CheckSeqDir ['y' or 'n']: yes or no for checking fwd or rev sequence
%      direction. If yes, then it will swap the sequence direction to the
%      fwd direction.
%
%  NOTE
%    Will fill in the CDR3 start and end fields of VDJ data with single or
%    multiple integer values. The findVDJmatch and findGeneMatch combo will
%    uses these seed locations to speed up the alignment process, and
%    updateVDJdata will finally select the final values.
%
%  See also findVDJmatch, findGeneMatch, updateVDJdata

function VDJdata = seedCDR3position(VDJdata,NewHeader,Xmap,X,Nleft,Nright,CheckSeqDir)
%Extract getHeaderVar variables since parfor can't handle it.
SeqLoc = findCell(NewHeader,{'nucleotide','Seq'});
CDR3startLoc = findCell(NewHeader,{'CDR3_Start'});
CDR3endLoc = findCell(NewHeader,{'CDR3_End'});

%Setup seed and inputs
Xseed = getGeneSeed(Xmap,X,Nleft,Nright,'nt');
CheckSeqDir = upper(CheckSeqDir(1));
X = upper(X);

%Setup the convolveSeq input structure for faster alignment
P.AllowedMiss = 0;
P.Alphabet = 'nt';
P.CheckSeq = 'yes';
P.DiagIdx = [];
P.ExactMatch = 'no';
P.PreferSide = 'none';
if strcmp(X,'V')
    P.TrimSide = 'right';   %will trim poorly matched regions in this side
    P.PreferSide = 'left';  %if tied alignments, will favor one towards left of seed and seq
    P.PenaltySide = 'left'; %unused nts of seed in this side will be unvaroable
    SpecSeed = 'TGT';       %conserved C
    CDR3colLoc = CDR3startLoc;  %Where in VDJdata to store the location
elseif strcmp(X,'J')
    P.TrimSide = 'left';    %will trim poorly matched regions in this side
    P.PreferSide = 'right'; %if tied alignments, will favor one towards left of seed and seq
    P.PenaltySide = 'right';%unused nts of seed in this side will be unvaroable
    SpecSeed = 'TGG';       %conserved W
    CDR3colLoc = CDR3endLoc;    %Where in VDJdata to store the location
end

%Find the CDR3start or CDR3end locations. Flip seq too if needed.
parfor j = 1:size(VDJdata,1)
    Tdata = VDJdata(j,:);  %Extract this to "slice" VDJdata
    Seq = Tdata{SeqLoc};
    
    %Do seed alignments, forward sense
    AlignScores = zeros(size(Xseed,1),3);
    for x = 1:size(Xseed,1)
       [Score,~,StartAt,MatchAt] = convolveSeq(Xseed{x},Seq,P);
       if StartAt(2) > 0
            AnchorLoc = Nleft - StartAt(2) + 2;
        else
            AnchorLoc = Nleft + abs(StartAt(2)) + 1;
       end       
       AlignScores(x,:) = [Score(1)/(diff(MatchAt)+1) Score(2) AnchorLoc];
    end

    %Do seed alignments on reverse complement sequence
    if CheckSeqDir == 'Y'
        AlignScoresR = zeros(size(Xseed,1),3);
        SeqNTR = seqrcomplement(Seq);
        for x = 1:size(Xseed,1)
            [Score,~,StartAt,MatchAt] = convolveSeq(SeqNTR,Xseed{x},P);
            AlignScoresR(x,:) = [Score(1)/(diff(MatchAt)+1) Score(2) StartAt(2)];
        end

        %Accept complement sequence if it's a higher score AND less positions
        if max(AlignScoresR(:,2)) > max(AlignScores(:,2))
            FwdPos = unique(AlignScores(:,3));
            RevPos = unique(AlignScores(:,3));
            if length(RevPos) <= length(FwdPos) %Means seeds are converging more on reverse direction          
                Tdata{1,SeqLoc} = SeqNTR;
                AlignScores = AlignScoresR;
            end
        end
    end
    
    %Finalize all potential CDR3start or CDR3end locations
    SpecPos = regexp(Seq,SpecSeed); %Special seed for V and J, which is the 'TGT' and the 'TGG'
    CDR3pos = unique([AlignScores(:,3); SpecPos(:)]);
    if X == 'J' 
        CDR3pos = CDR3pos + 2; %Need to include 2 nt of codon.
    end
    InvalidPos = (CDR3pos > length(Seq))  |  (CDR3pos < 1);
    CDR3pos(InvalidPos) = [];

    %Update to VDJdata ONLY if there is a seed
    if ~isempty(CDR3pos)
        Tdata{1,CDR3colLoc} = CDR3pos'; %Update Tdata instead of VDJdata b/c of parfor sliced-variable rule.
        VDJdata(j,:) = Tdata; %Add Tdata to sliced variable VDJdata
    end
end