%unique2 will find unique text and number in a cell array. This
%overcomes the issue when cells have mixed text and numbers, which cannot
%use MATLAB built-in 'unique.m' function. This is generally faster if
%dealing with finding indices (IC) for > 1000 sequences.
%
%  [U, IA, IB, IC] = unique2(C)
%
%  INPUT
%    C: Mx1 cell of text and numbers
%  
%  OUTPUT
%    U: unique cells of text and numbers
%    IA: index of A that will yield C,    C = A(IA)
%    IB: index of C for each entity of A, A = C(IB)
%    IC: cell of length(U) contain indices for each cell, U{1} = indices of
%      A for 1st unique cell in U
%
%  EXAMPLE
%    C = {'a' 'd' 'b' 'c' 'c' 'b'};
%    [U, IA, IB, IC] = unique2(C)
%    U =
%         'a'
%         'b'
%         'c'
%         'd'
%    IA =
%          1
%          3
%          4
%          2
%    IB =
%          1
%          4
%          2
%          3
%          3
%          2
%    IC =
%         [1]
%         [3;6]
%         [4;5]
%         [2]

function [Unq, IA, IB, IC] = unique2(C)
if isempty(C)
    Unq = C;
    IA = [];
    IB = [];
    IC = {};
    return
end

if isnumeric(C)
    C = num2cell(C);
end
[U, SortIdx] = sort2(C(:));

UnqIdx  = zeros(length(C), 1, 'uint64'); %Used to identify unique value
TempIdx = zeros(length(C), 1, 'uint64'); %Used for storing indices
IdxCell = cell(length(C), 1); %Used for storing cell of indices
CurVal = C{1};
UnqIdx(1) = 1;
t = 1;
c = 1;
u = 2;
for k = 1:length(U)
    if isequal(U{k}, CurVal)
        TempIdx(t) = k; 
        t = t+1;
        continue
    end
    IdxCell{c} = TempIdx(1:t-1);
    TempIdx(1) = k;
    t = 2;
    c = c + 1;
    CurVal = U{k};
    UnqIdx(u) = k;
    u = u + 1;
end
if all(U{k} == CurVal)
    IdxCell{c} = TempIdx(1:t-1);
end

Unq = U(UnqIdx(1:u-1));
IC = cellfun(@(x) SortIdx(x), IdxCell(1:c), 'unif', false);
IB = repelem([1:length(Unq)]', cellfun(@length, IC));
IB(SortIdx) = IB;
IA = SortIdx(UnqIdx(1:u-1));
