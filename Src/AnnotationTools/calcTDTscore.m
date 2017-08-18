%calcTDTscore will calculate the probability terminal deoxynucleotidyl
%transferase (TDT) assembles a nucleotide sequence, relative to random nt
%insertion. Returns output as ProbTDT/ProbRandom. It will check forwards
%and complement sequence to determine the score, which is the maximum score
%of the 2 direction.
%
%  TDTscore = calcTDTscore(Seq)
%
%  TDTscore = calcTDTscore(Seq,Prob)
%
%  INPUT
%    Seq: nucleotide sequence consisting of ACGTU letters
%    Prob: 1x4 matrix of probability that TDT adds a nucleotide. Must add
%      to 1, and will replace the default probability values.
%      Default values based on BRILIA analysis of mouse C57BL6 repertoire
%      Prob = [0.25 0.08 0.60 0.07]; %[Pa Pc Pg Pt] 
%
%  OUPUT
%    TDTscore: Calculated as P_TDT/(P_Random + P_TDT), where 
%      P_TDT = prod(Px,i) for all ith letter in sequence of length L
%      P_Random = 0.25^L
%
%  NOTES
%    - Both the forward and complement sequence are used to determine the
%      maximum TDTscore, which is the one that is returned. 
%    - If no sequence are placed inside, will return an empty score [].
%    - Ambiguous characters (eg 'N' or 'X') are NOT counted in the score.
%
%  See also findBetterD, trimGeneEdge

function TDTscore = calcTDTscore(SeqF,varargin)
%Set the probability matrix
if length(varargin) >= 1
    Prob = varargin{1}; 
    if length(Prob) ~= 4
        error('Prob should be a double value of length = 4');
    elseif max(Prob) > 1 || min(Prob) < 0
        error('Prob individual value must range from 0 to 1');
    elseif sum(Prob) ~= 1
        error('Sum of Prob should be equal to 1');
    end
else
    Prob = [0.25 0.08 0.60 0.07]; %[Pa Pc Pg Pt]
end

%If no sequence, return empty value
if isempty(SeqF)
    TDTscore = [];
    return
end

%Remove ambiguous characters and get +/- sense sequence
SeqFint = nt2int(SeqF);
DelThis = SeqFint > 4 | SeqFint == 0;
SeqF(DelThis) = [];
SeqFint(DelThis) = [];
SeqR = seqcomplement(SeqF);
SeqRint = nt2int(SeqR);

%Calculate the forward and reverse TDTscores
ProbF = prod(Prob(SeqFint));
ProbR = prod(Prob(SeqRint));
ScoreF = ProbF/(0.25^length(SeqF) + ProbF);
ScoreR = ProbR/(0.25^length(SeqF) + ProbR);

%Get only the larger score
TDTscore = max([ScoreF ScoreR]);
