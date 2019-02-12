%conformDist will conform multiple distribution data to have the same bin
%categories. This is used mainly to combine multiple frequency Mx2 cells
%into one Mx(N+1) cell, where the 1st column is the unique labels.
%
%  C = conformDist(D1, D2, ..., DN)
%
%  INPUT
%    Dn: Mn x 2 cell of label (col1) and frequency (col2) for nth dist.
%  
%  OUTPUT
%    C: M x (N+1) cell of frequencies, where C(:, 1) are the unique, sorted
%      labels from all distributions.
%
%  EXAMPLE
%    D{1} = {'A1' 10; 'A2' 29; 'A3' 40};
%    D{2} = {'A1' 20; 'A3' 40; 'A5' 60};
%    D{3} = {'A2' 30; 'A4' 30; 'A5' 4};
%    C = conformDist(D{:}); 
%    C = 
%         'A1'    [10]    [20]    [ 0]
%         'A2'    [29]    [ 0]    [30]
%         'A3'    [40]    [40]    [ 0]
%         'A4'    [ 0]    [ 0]    [30]
%         'A5'    [ 0]    [60]    [ 4]
%
function C = conformDist(varargin)
if any(~cellfun(@iscell, varargin))
    error('%s: Inputs must be MxN cells.', mfilename)
end

NumCol = sum(cellfun(@(x) size(x, 2)-1, varargin));

List = cellfun(@(x) x(:, 1), varargin, 'unif', false);
Label = sort2(unique2(vertcat(List{:})));
C = [Label(:) num2cell(zeros(length(Label), NumCol))];
k = 2;
for f = 1:length(varargin)
    if all(cellfun(@isnumeric, Label)) && all(cellfun(@isnumeric, List{f}))
        [~, IdxA, IdxB] = intersect(cell2mat(Label), cell2mat(List{f}), 'stable');
    else
        [~, IdxA, IdxB] = intersect(Label, List{f}, 'stable');
    end
    for j = 2:size(varargin{f}, 2)
        C(IdxA, k) = varargin{f}(IdxB, j);
        k = k + 1;
    end
end
