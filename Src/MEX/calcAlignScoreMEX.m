%calcAlignScoreMEX takes a logical MxN matrix returned from a sequence
%comparison resuls (SeqA == SeqB) and calculates the alignment score per
%each row using the following equation:
%SUM((consecutive matches)^2) - SUM((consecutive mismatches)^2) 
%
%  Score = calcAlignScoreMEX(MatchResults)
%
%  Score = calcAlignScoreMEX(MatchResults, AllowedMiss)
%
%  Score = calcAlignScoreMEX(MatchResults, AllowedMiss, PenaltySide)
% 
%  INPUT
%    MatchResults: 1xN logical array
%    AllowedMiss: number of single-0 gap (eg, 1 0 1) that can be ignored 
%      from left to right when calculating the consecutive match segment.
%      EX: 1 0 1, AllowedMiss = 1 will return Score = 2^2 = 4.
%      EX: 1 0 1, AllowedMiss = 0 will return Score = 1 - 1 + 1 = 1;
%    PenaltySide: 'n', 'r', 'l', 'b' subtract penalty^2 for mismatched edges
%
%  OUTPUT
%    Score: the alignment score
%
%  EXAMPLE
%    MatchResults = [1 0 1 0 1 0 1 1 1 1]>0;
%    AllowedMiss = 3;
%    Score = calcAlignScoreMEX(MatchResults, AllowedMiss);
%          =  49;
%    Score = calcAlignScoreMEX(MatchResults, 0);
%          =  16;
% 
%
