%unique2 is a wrapper for matlab's unique, which returns the cell array of
%indices of unique elements or rows, which saves times if the indices are
%needed.
%
%  [U, IU, IC, IdxCell] = unique2(C)
%
%  [U, IU, IC, IdxCell] = unique2(C, 'rows')
%
%  [U, IU, IC, IdxCell] = unique2(C, 'rows', varargin)
%
%  INPUT
%    C: matrix or cell array
%    'rows': sort C by rows
%    varargin: any other inputs accepted by unique.m
%  
%  OUTPUT
%    U: unique array of C
%    IU: index of U that will yield C, Ex: C = A(IU, :)
%    IC: index of C for each row of U, Ex: A = C(IC, :)
%    IdxCell: cell array of size(U,1)x1 containing indices for each unique 
%      row of U that is found in C.    Ex: U(1, :) = C(IdxCell{1}(1), :)
%
%  EXAMPLE
%    C = {'a' 'd' 'b' 'c' 'c' 'b'};
%    [U, IU, IC, IdxCell] = unique2_New(C)
%    U =
%         'a'
%         'b'
%         'c'
%         'd'
%    IU =
%          1
%          3
%          4
%          2
%    IC =
%          1
%          4
%          2
%          3
%          3
%          2
%    IdxCell =
%         [1]
%         [3;6]
%         [4;5]
%         [2]

function [Unq, IU, IC, IdxCell] = unique2(C, varargin)
if isempty(C)
    Unq = C;
    IU = [];
    IC = [];
    IdxCell = {};
    return
end

%Different cases to handle
RowOptLoc = strcmpi(varargin, 'rows');
UseSortRows = any(RowOptLoc);
if UseSortRows
    varargin = varargin(~RowOptLoc);
    if isnumeric(C)
        [U, SortIdx] = sortrows(C, varargin{~RowOptLoc});
    elseif iscell(C)
        [U, SortIdx] = sort2(C, varargin{~RowOptLoc}, 'rows');
    else
        error('%s: Cannot find unique rows of this class of array', mfilename);
    end
else
    C = C(:);
    [U, SortIdx] = sort2(C, varargin{:});
end

%Determine the indices
UnqIdx = repelem(0,  size(C, 1), 1);    %Used to identify unique value
TmpIdx = repelem(0,  size(C, 1), 1);    %Used for storing indices to prevent growing matrices
IdxCell = repelem({[]}, size(U, 1), 1); %Used for storing cell of indices
UnqIdx(1) = 1;
TmpIdx(1) = 1;
CurVal = U(1, :);
u = 1; %Counter for UnqIdx
t = 1; %Counter for TempIdx
i = 0; %Counter for IdxCell
for k = 2:size(U, 1)
    if isequal(U(k, :), CurVal)
        t = t + 1;
        TmpIdx(t) = k; 
    else
        i = i + 1;
        u = u + 1;
        IdxCell{i} = SortIdx(TmpIdx(1:t)); %Remember to get the original index of C
        CurVal = U(k, :);
        UnqIdx(u) = k;
        TmpIdx(1) = k; %Reset
        t = 1;         %Reset
    end
end
IdxCell{i+1} = SortIdx(TmpIdx(1:t)); %Save the last unique group, and remember to get the original index of C

Unq = U(UnqIdx(1:u), :);
IC = repelem((1:u)', cellfun('length', IdxCell(1:u)));
IC(SortIdx(:, 1)) = IC;
IU = SortIdx(UnqIdx(1:u));
IdxCell = IdxCell(1:u);