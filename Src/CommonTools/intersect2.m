%intersect2 returns cell array of indices of intersecting values, which
%differs from intersect that returns only the first occurence that
%intersects.
%
%  A = {'a' 'b' 'c' 'a' 'a'};
%  B = {'a' 'c'};
%  [C, IA, IB] = intersect2(A, B)
%  C =
%     'a'
%     'c'
%  IA =
%     [1; 4; 5]
%     [      3]
%  IB =
%     [1]
%     [2]

function varargout = intersect2(A, B, varargin)
[UnqA, ~, ~, IdxA] = unique2(A);
[UnqB, ~, ~, IdxB] = unique2(B);
[UnqC, IA, IB] = intersect(UnqA, UnqB);
varargout{1} = UnqC;
varargout{2} = IdxA(IA);
varargout{3} = IdxB(IB);