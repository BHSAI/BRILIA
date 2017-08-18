%calcSeqMutability will calculate the mutability of each nt based on
%published mutability index values.
%
%  MutVal = calcSeqMutability(Seq, mutMap)
%
%  INPUT
%    Seq: a sequences string to analyze. Ambiguous letters void
%      calculations at that position.
%    mutMap: hash map relating a tri-nucleotide string to its 1x3
%      mutability index for the 3 nt positoins. (getMutabilityMap.m)
%
%  OUTPUT
%    MutVal: a 1xN matrix of maximum mutability indices per each Seq nt
%
%  See also getMutabilityMap

function MutVal = calcSeqMutability(Seq, mutMap)
MutVal = zeros(3,length(Seq));
k = 1; %Shift the row to 1, 2, 3, 1, 2, ...
for j = 1:length(Seq)-2
    TriNuc = Seq(j:j+2);
    try
        MutVal(k, j:j+2) = mutMap(TriNuc);
        k = k + 1;
        if k > 3; k = 1; end
    catch
        fprintf('%s: Could not find mutability for %s.\n', mfilename, TriNuc);
    end
end
MutVal = max(MutVal, [], 1);
