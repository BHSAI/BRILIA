%calcAlignScore take a binary matrix of 1 or 0 and calculate the alignment
%score. This score higher for consecutive matches, versus identity matches.
%
%  Score = calcAlignScore(MatchResults), where MatchResults can either by a
%  binary matrix of match (1) and mismatch (0), OR, MatchResults can be a
%  3xN char matrix from an alignment results.
%
%  Score = calcAlignScore(MatchResults,AllowedMiss) will tolerate certain
%  mismatched pairs, from left to right, but only for point mutations.
%
%  The formula is SUM( (consecutive matches)^2 ) + SUM( (consecutive mismatches)^2 ) 
%
%  EX (2 point mutations)
%    Seq1 = 'GGGGGG'
%    Seq2 = 'GTGGTG'
%    MatchResults = Seq1 == Seq2;
%    Score1 = calcAlignScore(MatchResults,0)
%        Score1 = 4
%    Score2 = calcAlignScore(MatchResults,1)
%        Score2 = 9
%    Score3 = calcAlignScore(MatchResults,2)
%        Score2 = 16
%
%  EX (1 consecutive double mutations, so doesn't count)
%    Seq1 = 'GGGGGG'
%    Seq2 = 'GGTTGG'
%    MatchResults = Seq1 == Seq2;
%    Score1 = calcAlignScore(MatchResults,0)
%        Score1 = 4
%    Score2 = calcAlignScore(MatchResults,2)
%        Score2 = 4


function Score = calcAlignScore(MatchResults,varargin)
%Determine AllowedMismatch
if ~isempty(varargin)
    AllowedMiss = varargin{1};
else
    AllowedMiss = 0;
end

%Determine the per bp match score, depending the the Scoring Mode.
if size(MatchResults,1) == 1
    PosMatchVal = double(MatchResults);
else
    PosMatchVal = double(MatchResults(2,:) == '|');
end

%if AllowedMiss > 0, need to input leneiency rule on mismatches
if AllowedMiss > 0
    PosMatchIdx = find(PosMatchVal == 1);
    DiffPosIdx = find(diff(PosMatchIdx) == 2);
    if ~isempty(DiffPosIdx)
        if length(DiffPosIdx) > AllowedMiss
            DiffPosIdx = DiffPosIdx(1:AllowedMiss);
        end
        PosMatchVal(PosMatchIdx(DiffPosIdx)+1) = 2;
    end
end
NegMatchVal = 1-PosMatchVal;
   
%Determine the match score
Score = 0;
ConsHits = 0;
ConsMiss = 0;
MissPenalty = 0;
for j = 1:length(PosMatchVal)
    if PosMatchVal(j) > 0
        if PosMatchVal(j) <= 1
            ConsHits = ConsHits + PosMatchVal(j);
        else
            ConsHits = ConsHits - MissPenalty;
            MissPenalty = MissPenalty + (1-MissPenalty)/2;  %Increase penalty per miss.
        end
    elseif PosMatchVal(j) == 0
        Score = Score + ConsHits.^2;
        ConsHits = 0;
    end

    if NegMatchVal(j) > 0
        ConsMiss = ConsMiss + NegMatchVal(j);
    elseif NegMatchVal(j) == 0
        Score = Score - ConsMiss.^2;
        ConsMiss = 0;
    end
end

%Do final Addition and Subtraction.
Score = Score - ConsMiss.^2 + ConsHits.^2;