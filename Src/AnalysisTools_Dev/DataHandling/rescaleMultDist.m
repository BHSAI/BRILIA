%scaleMultDist will take the Mx(N+1) data cell and a grouping data
%(optional) to rescale the distribution with respect to the calculated mean
%of each data bin (along row). 
%
%  ScaleDist = scaleMultDist(C)
%
%  [ScaleDist, PerErr] = scaleMultDist(C, Group)
%
%  [ScaleDist, PerErr] = scaleMultDist(C, Group, RefGroup)
%
%  INPUT
%    C: conformed Mx(N+1) cell where 1st col is bin labels (conformDist.m)
%    Group: 1xN or Nx1 vector of integers that assigns a distribution to
%      a treatment group 
%    RefGroup: the group number to rescale all other bars to this
%      RefGroup's mean values
%
%  OUTPUT
%    ScaleDist: conformed and rescaled Mx(N+1) cell
%    PerErr: percentage error around group's mean
%
%  EXAMPLE
%    C = {'a' 10 20 50; 
%         'b' 20 40 100;
%         'c' 10 10 50};
%    Group = [1 1 2];
%    [ScaleDist, PerErr] = scaleMultDist(C, Group, 1)
%    ScaleDist =
%         'a'    [15.0003]    [14.9996]    [15.0000]
%         'b'    [30.0006]    [29.9992]    [30.0000]
%         'c'    [15.0003]    [ 7.4998]    [15.0000]
%    [ScaleDist, PerErr] = scaleMultDist(C, Group, 2)
%    ScaleDist =
%         'a'    [ 50.0010]    [49.9999]    [ 50.0050]
%         'b'    [100.0019]    [99.9999]    [100.0099]
%         'c'    [ 50.0010]    [25.0000]    [ 50.0050]
%
function [ScaleDist, PerErr, Scalor] = rescaleMultDist(C, Group, varargin)

if nargin == 1 
    Group = ones(1, size(C, 2)-1);
else
    if (size(C, 2)-1 ~= numel(Group))
        error('%s: size(C, 2)-1 must equal length(Group)', mfilename);
    end
end

UnqGrouping = unique(Group);
Freq = cell2mat(C(:, 2:end));
ScaleFreq = zeros(size(Freq));
ErrVal  = zeros(size(Freq));
for k = 1:length(UnqGrouping)
    Idx = find(UnqGrouping(k) == Group);
    MeanFreq = mean(Freq(:, Idx), 2);
    for j = 1:length(Idx)
        ScaleFreq(:, Idx(j)) = rescaleDist(MeanFreq, Freq(:, Idx(j)));
        ErrVal(:, Idx(j)) = ( ScaleFreq(:, Idx(j)) - MeanFreq ) ./ MeanFreq;
    end
end

ScaleDist = C;
ScaleDist(:, 2:end) = num2cell(ScaleFreq);

if nargout >= 2
    PerErr = C;
    PerErr(:, 2:end) = num2cell(ErrVal);
end

Scalor = zeros(1, length(Group));
if ~isempty(varargin)
    ScaleDist = C;
    RefIdx = (Group == varargin{1});
    RefDist = mean(Freq(:, RefIdx), 2);
    for j = 1:size(Freq, 2)
        [SDist, Scalor(j)] = rescaleDist(RefDist, Freq(:, j));
        ScaleDist(:, j+1) = num2cell(SDist); 
    end
end
