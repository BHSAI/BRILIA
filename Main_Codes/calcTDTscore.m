%calcTDTscore will take a nucleotide sequence, Seq, and TDTmatrix to
%compute the likelihood the sequence originates from TDT nucleotide
%insertion versis random mutation. Returns output as ProbTDT/ProbRandom.

function TDTscore = calcTDTscore(SeqF) %,TDTadjMatrix,TDTjoinMatrix)
Prob = [0.25 0.08 0.60 0.07];

if isempty(SeqF)
    TDTscore = 1;
    return
end

SeqFint = nt2int(SeqF);
DelThis = SeqFint > 4;
SeqF(DelThis) = [];
SeqFint(DelThis) = [];

SeqR = seqcomplement(SeqF);
SeqRint = nt2int(SeqR);
SeqRint(SeqRint>4) = [];

ProbF = 1; 
ProbR = 1;

for j= 1:length(SeqF)
    ProbF = ProbF*Prob(SeqFint(j));
    ProbR = ProbR*Prob(SeqRint(j));
end

ScoreF = ProbF/(0.25^length(SeqF) + ProbF);
ScoreR = ProbR/(0.25^length(SeqF) + ProbR);

TDTscore = max([ScoreF ScoreR]);