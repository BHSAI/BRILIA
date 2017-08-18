%calcAlignScore takes a logical MxN matrix returned from a sequence
%comparison resuls (SeqA == SeqB) and calculates the alignment score per
%each row using the following equation:
%SUM( (consecutive matches)^2 ) + SUM( (consecutive mismatches)^2 )
%
%  Score = calcAlignScore(MatchResults)
%
%  Score = calcAlignScore(MatchResults,AllowedMiss)
%
%  Score = calcAlignScore(MatchResults,AllowedMiss,GoodIdx)
%
%  INPUT
%    MatchResults: MxN logical matrix of match (1) and mismatch (0), OR a
%      3xN char matrix with 2nd row having "|" as matches
%    AllowedMiss: Integer number of allowed point misses allowed, from left
%      to right, in order to elongate a consecutive match segment. Score
%      for [1 0 1] is 1 with AllowedMiss = 0, and 4 if AllowedMiss = 1;
%    GoodIdx: Empty or MxN logical matrix indiciating which nt matches
%      should be counted towards the score. Used mainly if you don't want
%      edges to be counted towards the scoring function. If GoodIdx is
%      empty, score will be calculated based on location of first and last
%      match.
%
%  OUTPUT
%    Score: Mx1 alignment scores
%
%  EXAMPLES
%    Case 1) 2 point mutations
%      Seq1 = 'GGGGGG'
%      Seq2 = 'GTGGTG'
%      MatchResults = Seq1 == Seq2;
%      Score1 = calcAlignScore(MatchResults,0)
%          Score1 = 4
%      Score2 = calcAlignScore(MatchResults,1)
%          Score2 = 9
%      Score3 = calcAlignScore(MatchResults,2)
%          Score2 = 16
%
%    Case 2) 1 consecutive double mutations (so must accept mismatch)
%      Seq1 = 'GGGGGG'
%      Seq2 = 'GGTTGG'
%      MatchResults = Seq1 == Seq2;
%      Score1 = calcAlignScore(MatchResults,0)
%          Score1 = 4
%      Score2 = calcAlignScore(MatchResults,2)
%          Score2 = 4
%
%    Case 3) Using GoodIdx input
%      SeqA = 'ABCD'
%      SeqB = 'CDFADSD';
%      [Diag1, Diag2, GoodIdx] = makeDiagonalSeq(SeqA,SeqB)
%      Scores = calcAlignScore(Diag1 == Diag2,GoodIdx)
%
%  See also alignSeq, makeDiagonalSeq

function AllScores = calcAlignScore(MatchResults,varargin)
%Determine if you have a binary matrix, or a 3xN alignment char.
if ischar(MatchResults)
    if size(MatchResults == 3)
        MatchResults = MatchResults(2,:) == '|';
    end
end
    
%Initlized or get the extra inputs
AllowedMiss = 0;
GoodIdx = [];
if length(varargin) >= 1
    AllowedMiss = varargin{1};
    if length(varargin) >= 2
        GoodIdx = varargin{2};
    end
end

%Create the matrix of action, 1 = match, 2 = allowed miss
ActionMat = int8(MatchResults); %Elongate and add to score^2
if AllowedMiss > 0 %Apply leniency rule, left to right
    MissLocs = MatchResults(:,1:end-2) & MatchResults(:,3:end) & ~MatchResults(:,2:end-1);
    CumSumMiss = cumsum(MissLocs,2);
    AllowedMissLocs = 2*int8(MissLocs & (CumSumMiss <= AllowedMiss));
    ActionMat(:,2:end-1) = ActionMat(:,2:end-1) + AllowedMissLocs; %Elongate but do not increase score
end

%Mark places NOT to calculate scores
if ~isempty(GoodIdx)
    ActionMat(~GoodIdx) = -1; %Do not count these as misses
    MinAction = 0; %Score counts from 1st to last occurence of GoodIdx
else
    MinAction = 1; %Score counts from 1st to last occurence of MatchResults
end

%Compute the scores now
AllScores = zeros(size(MatchResults,1),1);
for k = 1:size(MatchResults,1)
    CurAction = ActionMat(k,ActionMat(k,:)>=0); %Remember to ignore bad idx, which have negative action value.    
    CurActionLoc = find(CurAction >= MinAction);

    %Possible to have no match
    if isempty(CurActionLoc)
        continue
    end
    StartLoc = CurActionLoc(1);
    EndLoc = CurActionLoc(end);
    
    Score = 0;
    ConsHits = 0;
    ConsMiss = 0;
    for j = StartLoc:EndLoc %StartLoc(1):size(CurAction,2)
        if CurAction(j) > 0;
            ConsHits = ConsHits + (CurAction(j) == 1);
            if ConsMiss > 0
                Score = Score - ConsMiss.^2;
                ConsMiss = 0;
            end
        else
            ConsMiss = ConsMiss + 1;
            if ConsHits > 0
                Score = Score + ConsHits.^2;
                ConsHits = 0;
            end
        end            
    end
    
    %Perform final additions to AllScores
    AllScores(k) = Score + ConsHits.^2 - ConsMiss.^2;
end
