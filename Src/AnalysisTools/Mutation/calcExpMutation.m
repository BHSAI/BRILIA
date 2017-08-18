%calcExpMutation will calculate the expected replacement AA mutation given
%a NT sequence.


function [Rexp, varargout] = calcExpMutation(Seq, mutMap, codonSRMap, relMutMat)
%Max mutability index per position
Fgl = calcSeqMutability(Seq, mutMap);

%Rel mutability per position i and from x to y
IntSeq = nt2int(Seq);
Mgl = zeros(4, length(Seq));
for j = 1:4
    Idx = IntSeq == j;
    Mgl(:, Idx) = repmat(relMutMat(j, :)', 1, sum(Idx));
end

%Binary map of if a -> b leads to silent mutation (0) or replacement (1)
Igl = zeros(4, length(Seq));
for j = 1:3:length(Seq)-2
    Igl(:,j:j+2) = codonSRMap(Seq(j:j+2));
end

%Compute the expected replacement and silent mutation
Rexp = sum(sum(repmat(Fgl, 4, 1) .* Mgl .* Igl));
if nargout == 2
    Sexp = sum(sum(repmat(Fgl, 4, 1) .* Mgl .* (1-Igl)));
    varargout{1} = Sexp;
end
