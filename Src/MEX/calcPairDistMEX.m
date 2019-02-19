%calcPairDistMEX computes the hamming distance and number of valid hotspot 
%mutations from SeqA to SeqB, and vice versa. A valid hotspot mutation is
%one that agrees with the motif (WRC, GYW, WA, TW) AND has a pairwise 
%substitution pattern that agrees with deaminase subsitution patterns. 
%A valid Motif is +1, and a valid pairwise mutation is +1, for a maximum 
%point of 2 for classical hotspot mutation. Worse is -2 for disagreeing 
%motif AND pairwise mutation. If [valid mut] + [invalid mut] < 0, then SeqA 
%cannot be parent of SeqB. It is possible SeqA->SeqB, and SeqB->SeqA return
%negative ShmDir scores, as in SeqA and SeqB are not linkable under the same
%clonotype. Consecutive mismatches add M^2 score, decrease chance SeqA and 
%SeqB are related.
%
%    N  A  C  G  T   Interpretted as:
%  ---------------   --------------
%N | 0  0  0  0  0   N -> - A C G T
%A | 0  0  0  1  1   A -> N - C G T
%C | 0 -1  0 -1  1   C -> N A - G T 
%G | 0  1 -1  0 -1   G -> N A C - T
%T | 0 -1  0 -1  0   T -> N A C G -
% 
%  [Ham, Motif, Mut, Penalty, ShmDist] = calcPairDistMEX(SeqA, SeqB);
% 
%  [Ham, Motif, Mut, Penalty, ShmDist] = calcPairDistMEX(SeqList);
%
%  INPUT
%    SeqA: nt sequence string
%    SeqB: nt sequence string
%    SeqList: Mx1 cell of nt sequence strings
%
%  OUTPUT
%    Ham: MxM hamming distance matrix between SeqA and SeqB
%    Motif: MxM total count of hotspot (+1) and non-hotspot motifs (-1)
%    Mut: MxM total count of valid (+1) and invalid (-1) nt subsitutions
%    Penalty: MxM of total of Mi^2, where Mi = ith mismatch segment
%    ShmDist: MxM BRILIA's SHM distance. 
%
%  NOTE
%    Each row is the parent, each column is the child.
% 
%    The "SHM Distance" is calculated as:
%      ShmDist = Ham - (Motif + Mut - Penalty)/4;
% 
%  EXAMPLE
%
%    Seq{1} = 'AACAACGAACGAA';
%    Seq{2} = 'AATAATTAATGAA';
%    [Ham, Motif, Mut, Penalty, ShmDist] = calcPairDistMEX(Seq);
%    Ham =
%         0     4
%         4     0
%    Motif =
%         0     2
%         2     0
%    Mut =
%         0     2
%        -1     0
%    Penalty =
%         0     4
%         4     0
%    ShmDist =
%         0     6.00
%      6.75     0.00
% 
%
