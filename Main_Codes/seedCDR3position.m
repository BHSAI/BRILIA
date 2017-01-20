%seedCDR3position will find potential locations where 104C's 1st or 118W's
%last nt codon could reside in a sequence. This could also check both
%forward and reverse direction, returning the direction with the highest
%seed alignment score.
%
%  VDJdata = seedCDR3position(VDJdata,NewHeader,Xmap,X,Nleft,Nright,CheckSeqDir)
%
%  INPUT
%    Xmap [Vmap or Dmap]: Database for V or J germline gene
%    X ['V' or 'J']: Specificy what database gene is used
%    Nleft: number of nts left of 104C to use as alignment seed.
%    Nright: number of nts right of 104C to use as alignment seed.
%    CheckSeqDir ['y' or 'n']: yes or no for checking fwd or rev sequence
%      direction. If yes, then it will swap the sequence direction to the
%      fwd direction.

function VDJdata = seedCDR3position(VDJdata,NewHeader,Xmap,X,Nleft,Nright,CheckSeqDir)
%Bring some getHeaderVar variables out, since parfor can't handle it.
SeqLoc = findCell(NewHeader,{'nucleotide','Seq'});
CDR3startLoc = findCell(NewHeader,{'CDR3_Start'});
CDR3endLoc = findCell(NewHeader,{'CDR3_End'});

%Setup the rapid convolveSeq input structure
P.AllowedMiss = 0;
P.Alphabet = 'nt';
P.CheckSeq = 'yes';
P.ExactMatch = 'no';
P.DiagIdx = [];
P.PreferSide = 'none';
P.PenaltySide = 'none';
if strcmpi(X,'V')
    P.TrimSide = 'right';
    P.PreferSide = 'left';
elseif strcmpi(X,'J')
    P.TrimSide = 'left';
    P.PreferSide = 'right';
end
CheckSeqDir = CheckSeqDir(1);

%Get the seed sequence
Xseed = getGeneSeed(Xmap,X,Nleft,Nright,'nt');

parfor j = 1:size(VDJdata,1)
    Tdata = VDJdata(j,:);
    SeqNT = Tdata{SeqLoc};
    
    %Do seed alignments, forward sense
    AlignScores = zeros(size(Xseed,1),3);
    for v = 1:size(Xseed,1)
       [Score,~,StartAt,MatchAt] = convolveSeq(Xseed{v},SeqNT,P);
       if StartAt(2) > 0
            AnchorLoc = Nleft - StartAt(2);
        else
            AnchorLoc = Nleft + abs(StartAt(2)) + 1;
        end
        
        AlignScores(v,:) = [Score(1)/(diff(MatchAt)+1) Score(2) AnchorLoc];
    end

    %Do seed alignments on reverse complement sequence
    if strcmpi(CheckSeqDir,'y')
        AlignScoresR = zeros(size(Xseed,1),3);
        SeqNTR = seqrcomplement(SeqNT);
        for v = 1:size(Xseed,1)
            [Score,~,StartAt,MatchAt] = convolveSeq(SeqNTR,Xseed{v},P);
            AlignScoresR(v,:) = [Score(1)/(diff(MatchAt)+1) Score(2) StartAt(2)];
        end

        %Accept complement sequence if it's a higher score
        if max(AlignScoresR(:,2)) > max(AlignScores(:,2))
            Tdata{1,SeqLoc} = SeqNTR;
            AlignScores = AlignScoresR;
        end
    end

%     %Find the best score that also has 80% match, or sequences that that
%     %have 95% match and alignment score > (75% of seed length)^2. Should be
%     %stringent as you want good seeds. If it fails, all that happens is a
%     %full search thats slightly slower.
%     BestScore = max(AlignScores(:,2));
%     BestLocs = ((AlignScores(:,2) == BestScore) & (AlignScores(:,1) > 0.80)) | (AlignScores(:,1) >= 0.95 & AlignScores(:,2) > (0.75*length(Xseed))^2);
%     ModeLocs = (AlignScores(:,3) == mode(AlignScores(:,3))) & AlignScores(:,2) > 0;
%     AllLocs = ModeLocs | BestLocs;
% 
%     if max(BestLocs) == 0; continue; end %Don't do anything.

    %Save ANY positions that matches to seed. findGeneMatch will pick
    %correct seed location later. If you try to select seed now, you will
    %lose key locations.
    if strcmpi(X,'V') %for V, fill in CDR3 start positions
        CDR3pos = unique(AlignScores(:,3));
        CDR3pos(CDR3pos > length(SeqNT)) = [];
        Tdata{1,CDR3startLoc} = CDR3pos';
    else %For J, fill in CDR3 end positions
        CDR3pos = unique(AlignScores(:,3) + 2); %Need to include 2 nt of codon.
        CDR3pos(CDR3pos > length(SeqNT)) = [];
        Tdata{1,CDR3endLoc} = CDR3pos';
    end
    VDJdata(j,:) = Tdata;
end