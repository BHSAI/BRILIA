%findGeneMatch will look for the best sequence match to a gene database
%sequence set. This is the core alignment-based gene search.
%
%  GeneMatch = findGeneMatch(Seq, Xmap, X, MissRate)
%
%  GeneMatch = findGeneMatch(Seq, Xmap, X, MissRate, CDR3anchor)
%
%  GeneMatch = findGeneMatch(Seq, Xmap, X, MissRate, CDR3anchor, 'ForceAnchor')
% 
%  INPUTS
%    Seq: Character sequence
%    Xmap: Database DB.[X]map, where [X] is V, D, J, Vk, Vl, Jk, or Jl.
%    X ['V', 'D', 'J']: character specifying what segment to match. Does not
%      have to specify heavy/light chain Vs or Js.
%    MissRate [Integer >= 0]: How many missed nt/bp are allowed for
%      elongating a contiguous C region for alignment score. Ex: MissRate
%      0.02 = 2% of shorter seq length.
%    MissRate: How many point mutations per nt are allowed during
%      alignment. Default is 0.00, and ranges from 0.00 to 1.00.
%    CDR3anchor: Positions of potential 104C 1st nt of codon or 118W 3rd nt
%      of codon. CDR3 anchor can be a matrix of values, and providing it
%      speeds up alignment and reduces chance of incorrect J alignments, 
%      especially when dealing with shorter sequences cut at the 118W.
%   'ForceAnchor': will force findGeneMatch to find a solution that matches
%      to at least one of the anchor points. Will not attempt a unseeded
%      alignment if the seeded alignment gives a low alignment % (below 70%
%      match).
%
%  OUTPUT
%    GeneMatch: 1x6 cell array of cells containing gene match data as follows
%      Col1   Gene number in Xmap
%      Col2   Full gene name(s)
%      Col3   [LeftLength MiddleLength RighLength] of the reference gene
%      Col4   [LeftLength MiddleLength RighLength] of the Seq
%      Col5   [(# of matches) AlignmentScore]
%      Col6   3xN character alignment results
%
%  See also alignSeq, findVDJmatch, findVJmatch

function GeneMatch = findGeneMatch(Seq, Xmap, X, MissRate, varargin)
%Make sure Seq is a cell.
if iscell(Seq)
    if length(Seq) > 1
        warning('findGeneMatch: Seq should be a character. Taking only the first sequence of the cell.');
    end
    Seq = Seq{1};
end

X = upper(X(1)); %Want to handle only uppercase X

%Determine if 'ForceAnchor' option was given
ForceAnchor = 0;
ForceAnchorLoc = findCell(varargin, 'ForceAnchor', 'MatchCase', 'any');
if ForceAnchorLoc > 0
    ForceAnchor = 1;
    varargin(ForceAnchorLoc) = [];
end

%Determine if a valid, nonzero CDR3anchor was given.
CDR3anchor = 0;
if ~isempty(varargin) && ~isempty(varargin{1})
    CDR3anchor = sort(varargin{1});
end

%Override ForceAnchor. Can't enforce an anchor that does not exist.
if max(CDR3anchor) == 0 && ForceAnchor == 1
    ForceAnchor = 0;
end

%Ensure miss rate is valid
if MissRate < 0; MissRate = 0; end
if MissRate > 1; MissRate = 1; end

%Setup the input structure for alignSeq, which is faster.
P.MissRate = MissRate;
P.Alphabet = 'nt';
P.CheckSeq = 'yes'; %Assumes ambig char taken char of before this.
P.DiagIdx = [];
if X == 'V'
    P.TrimSide = 'right';
    P.PreferSide = 'left';
    P.PenaltySide = 'left';
elseif X == 'J'
    P.TrimSide = 'left';
    P.PreferSide = 'right';
    P.PenaltySide = 'right';
elseif X == 'D'
    P.TrimSide = 'both';
    P.PreferSide = 'none';
    P.PenaltySide = 'none';
else
    P.TrimSide = 'none';
    P.PreferSide = 'none';
    P.PenaltySide = 'none';
end
if max(CDR3anchor) > 0
    P.ExactMatch = 'yes';
else
    P.ExactMatch = 'no';
end    

%--------------------------------------------------------------------------
%Find the best ref gene map.

%If CDR3 anchors are provided, try to use this first to find gene match.
if P.ExactMatch(1) == 'y'
    %Determine best fit ref gene number via NT matching
    FitScore1 = -Inf*ones(length(CDR3anchor), size(Xmap, 1)); %Match/mismatch
    FitScore2 = -Inf*ones(length(CDR3anchor), size(Xmap, 1)); %Align Score
    InvalidMap = zeros(1, size(Xmap, 1), 'logical');
    for x = 1:size(Xmap, 1)
        %Determine if Xmap entry and anchor are viable
        if isempty(Xmap{x, 1}); 
            InvalidMap(x) = 1;
            continue; 
        end %Deleted reference seq       
        if Xmap{x, end} == 0; 
            InvalidMap(x) = 1;
            continue; 
        end %No valid anchor available
        
        if X == 'V'
            Xanchor = length(Xmap{x, 1}) - Xmap{x, end} + 1; %CDR3start location (1st nt of codon)
        else 
            Xanchor = Xmap{x, end} + 2; %CDR3end location (3rd nt of codon)
        end
        
        %Determine alignment score for Xmap and anchor points
        for q = 1:length(CDR3anchor)
            [SeqA, SeqB] = padtrimSeq(Seq, Xmap{x, 1}, CDR3anchor(q), Xanchor, 'min', 'min');
            [ScoreT, ~, StartAt, MatchAt] = alignSeq(SeqA, SeqB, P); 
            %Adjust Score1 based on where SeqA start/ends for V/J gene
            if StartAt(2) < 0
                SeqAstart = abs(StartAt(2)) + 1;
                SeqAend = SeqAstart + length(SeqA) - 1;
            else
                SeqAstart = 1;
                SeqAend = SeqAstart + length(SeqA) - 1;
            end
            if X == 'V'
                if MatchAt(1) > SeqAstart %SeqA left side was not matched. Penalize.
                    MatchAt(1) = SeqAstart;
                end
            elseif X == 'J'
                if MatchAt(2) < SeqAend %SeqA right side was not matched. Penalize.
                    MatchAt(2) = SeqAend;
                end
            end
      
            %Save the alignment scores
            FitScore1(q, x) = ScoreT(1)/(diff(MatchAt)+1);
            FitScore2(q, x) = ScoreT(2);                                
        end
    end
    
    %Find the best anchor point and Xmap number
    [BestRow, BestXmapNum] = find(FitScore2 == max(max(FitScore2(:, ~InvalidMap)))); %The index number in Xmap of best matches
    BestAnchor = CDR3anchor(BestRow);
    if length(BestRow) > 1 %If multiple anchors work, then...
        if X == 'V' %Select the furthest right start loc
            BestXmapNum = BestXmapNum(end, :); 
            BestAnchor = CDR3anchor(BestRow(end));
            BestRow = BestRow(end);
        elseif X == 'J' %Select the furthest right start loc
            BestXmapNum = BestXmapNum(1, :); 
            BestAnchor = CDR3anchor(BestRow(1));
            BestRow = BestRow(1);
        end
    end

    %In case you have to abandon a bad seed, just do full match
    BestIdentity = max(FitScore1(BestRow, BestXmapNum));
    if (BestIdentity < 0.70 && ForceAnchor == 0) || BestIdentity == 0
        P.ExactMatch = 'no';
    end
end
    
%If no anchor point, or previous one gave bad results, do full alighmnet.
if P.ExactMatch(1) == 'n'
    %Generate the DiagMatrixIdx here that is used by makeDiagonalSeq in
    %convolveSeq. This will speedup the gene matching process.
    
    %Determine the number of rows in the diagonal matrix
    SeqXmaxLen = 0;
    for x = 1:size(Xmap, 1)
        SeqXlen = length(Xmap{x, 1});
        if SeqXlen > SeqXmaxLen
            SeqXmaxLen = SeqXlen;
        end
    end
    Tlen = length(Seq)+SeqXmaxLen-1; %Length of the untrimmed diag matrix
    P.DiagIdx = repmat([-length(Seq)+2:1], Tlen, 1) + repmat([0:Tlen-1]', 1, length(Seq)); %Determine the index for Seq2  %SLOWEST    

    %Determine best fit ref gene number via NT matching
    FitScore = -Inf*ones(1, size(Xmap, 1));
    for x = 1:size(Xmap, 1)
        if isempty(Xmap{x, 1}); continue; end
        ScoreT = alignSeq(Seq, Xmap{x, 1}, P);
        FitScore(1, x) = ScoreT(2);
    end
    BestXmapNum = find(FitScore == max(FitScore)); %The index number in Xmap of best matches
    BestAnchor = 0;
end
%--------------------------------------------------------------------------

%Assemble the match cell, to rank by lowest ref gene deletion count
AlignResults = cell(length(BestXmapNum), 6);
for k = 1:length(BestXmapNum)
    w = BestXmapNum(k);
    SeqA = Seq;
    SeqB = Xmap{w, 1};
    Aadj = [0 0]; %Will always be negative to indicate how many nts were removed from left or right side of Seq A.
    Badj = [0 0]; %Will always be negative to indicate how many nts were removed from left or right side of Seq A.
    if P.ExactMatch(1) == 'y'
        if X == 'V'
            Xanchor = length(Xmap{w, 1}) - Xmap{w, end} + 1;
        else
            Xanchor = Xmap{w, end} + 2;
        end
        [SeqA, SeqB, Aadj, Badj] = padtrimSeq(SeqA, SeqB, BestAnchor, Xanchor, 'min', 'min'); %Using min, min option is key to ensuring negative Aadj and Badj values!
    end
    [AlignScore, Alignment, StartAt, MatchAt] = alignSeq(SeqA, SeqB, P); 

    %Calculate Left, Middle, and Right segment lengths. LMRsamp and LMRgerm
    if StartAt(2) < 0 %Seq2 starts first
        Ls = MatchAt(1) - abs(StartAt(2)) - 1 - Aadj(1);
        Ms = MatchAt(2) - MatchAt(1) + 1;
        Rs = abs(StartAt(2)) + length(SeqA) - MatchAt(2) - Aadj(2);

        Lg = MatchAt(1) - 1 - Badj(1);
        Mg = Ms;
        Rg = length(SeqB) - MatchAt(2) - Badj(2);
    else %Seq1 starts first
        Ls = MatchAt(1) - 1 - Aadj(1);
        Ms = MatchAt(2) - MatchAt(1) + 1;
        Rs = length(SeqA) - MatchAt(2) - Aadj(2);
        
        Lg = MatchAt(1) - StartAt(2) - Badj(1);
        Mg = Ms;
        Rg = StartAt(2) + length(SeqB) -1 - MatchAt(2) - Badj(2);
    end

    %Fill in AlignResults
    AlignResults{k, 1} = w; %GeneNum in the Xmap
    AlignResults{k, 2} = Xmap{w, 2}; %GeneName
    AlignResults{k, 3} = [Lg Mg Rg]; %Ref seq left del, middle length, right del
    AlignResults{k, 4} = [Ls Ms Rs]; %Sample seq left del, middle length, right del
    AlignResults{k, 5} = AlignScore(:)'; %Alignment score [SeqIdentity ConsecAlignScore]
    AlignResults{k, 6} = Alignment; %3xN cell matrix alignment results
end
AlignResults = condenseGeneMatch(AlignResults);

%Now find the match with the MINIMUM total deletion in the deletable ends
LMRgene = cell2mat(AlignResults(:, 3));
if X == 'V'
    MinLoc = find(LMRgene(:, 3)==min(LMRgene(:, 3)));
elseif X == 'J'
    MinLoc = find(LMRgene(:, 1)==min(LMRgene(:, 1)));
else
    Max5and3 = max(LMRgene(:, [1 3]), [], 2);
    MinLoc = find(Max5and3 == min(Max5and3));
end

%Report only 1st results if there are still ties.
GeneMatch = AlignResults(MinLoc(1), :);

%Make sure it's a valid match, else report null
if min(GeneMatch{4}) < 0 %Cannot have a negative segment. Caused by no match.
    GeneMatch = cell(1, 6);
end
