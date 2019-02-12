%getKmer will split a sequence in K-mers and store them as cells. This is
%used mainly before another function to count N-mer occurrences in a
%sequences.
%
%  Kmer = getKmer(Seq, N)
%
%  Kmer = getKmer(Seq, N, 'unique')
%
%  INPUT
%    Seq: cell array of sequences
%    N: N-mer to extract
%    'unique': will condense K-mers to unique K-mers
% 
%  OUTPUT
%    Kmer: cell array of cells of Kmers
%
%  EXAMPLE
%    Seq = {'ACGACGACG', 'AGGTGTAGTG'};
%    Kmer = getKmer(Seq, 3);
%    UnqKmer = getKmer(Seq, 3, 'unique')

function Kmer = getKmer(Seq, K, Option)
GetUniqueOnly = nargin == 3 && strcmpi(Option, 'unique');

Kmer = cell(size(Seq));
parfor j = 1:numel(Seq)
    MaxN = length(Seq{j}) - K + 1;
    Kmer{j} = cell(MaxN, 1);
    for k = 1:MaxN
        Kmer{j}{k, 1} = Seq{j}(k:k+K-1);
    end
    if GetUniqueOnly
        Kmer{j} = unique(Kmer{j});
    end
end
