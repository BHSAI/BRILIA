%countNT will count how many ACGTN nucleotide bases there are in a cell
%array of sequences of the same lengths.
%
%  Count = countNT(SEq)
%
%  INPUT
%    Seq: cell or char array of DNA/RNA sequences
%
%  OUTPUT
%    Count: 5xN matrix of number of base pairs observed. Rows 1-5 are
%    repsectively A, C, G, T, N.
%
%  EXAMPLE
%    Seq{1} = 'ACGTGANNCGTGTTT';
%    Seq{2} = 'ACGTTACACGACTTC';
%    Count = countNT(Seq);
%    Count =
%      2  0  0  0  0  2  0  1  0  0  1  0  0  0  0   %A
%      0  2  0  0  0  0  1  0  2  0  0  1  0  0  1   %C
%      0  0  2  0  1  0  0  0  0  2  0  1  0  0  0   %G
%      0  0  0  2  1  0  0  0  0  0  1  0  2  2  1   %T
%      0  0  0  0  0  0  1  1  0  0  0  0  0  0  0   %N or all others
%
function Count = countNT(Seq)

if iscell(Seq)
    Seq = vertcat(Seq{:});
elseif ~ischar(Seq)
    error('%s: Input must be a char or cell array of equal-length sequences', mfilename);
end
Seq = upper(Seq);

Count = zeros(5, size(Seq, 2));
Letters = 'ACGT';
for j = 1:length(Letters)
    Count(j, :) = sum(Seq == Letters(j), 1);
end
Count(5, :) = size(Seq, 1) - sum(Count(1:4, :), 1);