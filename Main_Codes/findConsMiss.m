%findConsMiss return the frequency count of a consensus mismatch of a
%seq or classifier array.
%
%  [ConsSeq, ConsCount] = findConsMiss(Seq,Type)
%  will take a trimmed and aligned Seq cell or char array, and return the
%  consensus mismatch count withing Seq and the absolute frequency of mismatch. Seq could be a NT
%  or Classifier sequence.
%    Mode = 'seq' will return consensus seq. If there is a tie, will return
%    an "n" in the nucleotide. Mode = 'classifier will return the matching
%    of classifiers instead, based on the upper case.

function [ConsSeq, ConsCount] = findConsMiss(Seq,Type)
%Convert to character array
if iscell(Seq)
    Seq = char(Seq);
end

ConsCount = zeros(1,size(Seq,2));
if strcmpi(Type,'seq')
    ConsSeq = seqconsensus(Seq);
    for k = 1:size(Seq,1)
        ConsCount = ConsCount + (Seq(k,:) ~= ConsSeq);
    end
elseif strcmpi(Type,'classifier')
    ConsSeq = ''; %Not applicable for classifiers
    for k = 1:size(Seq,1)
        ConsCount = ConsCount + isstrprop(Seq(k,:),'lower');
    end
end