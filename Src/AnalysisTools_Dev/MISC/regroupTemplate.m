%regroupTemplate will sum up template according to either a Nx1 cell or
%matrix of groupable contents, or a Mx1 cell of indices of Template matrix
%to group. 
%
%  GroupedTemplate = regroupTemplate(Template, Group)
%
%  [GroupedTemplate, Idx] = regroupTemplate(Template, Idx)
%
%  INPUT
%    Template: Nx1 matrix of numeric values
%    Group: Nx1 cell of values or char, where unique values are used to 
%      determine how to sum Template
%    Idx: Mx1 cell of Q-element indices to group up Template, normally from 
%      the 4th output of unique2.m   
%
%  OUTPUT
%    GroupedTemplate: Re-summed Template according to Group or Idx
%    Idx: cell of indices of the group, returned by unique2 4th output
%
%  EXAMPLE
%    Template = [1:10]';
%    Group = repmat({'a' 'a'; 'a' 'b'}, 5, 1)
%    [GroupedTemplate, Idx] = regroupTemplate(Template, Group)
%    Idx =
%      [5×1 double]  %1:2:9
%      [5×1 double]  %2:2:10

function [T, Idx] = regroupTemplate(Template, GroupOrIdx)
if iscell(GroupOrIdx) && all(cellfun('isclass', GroupOrIdx(:), 'double')) && numel(GroupOrIdx) <= numel(Template) %The Mx1 cell of grouped index was given. do no use unique2.
    Idx = GroupOrIdx;
else
    [~,~,~,Idx] = unique2(GroupOrIdx, 'rows');
end
T = zeros(numel(Idx), 1);
if iscell(Template)
    for j = 1:numel(Idx)
        T(j) = sum([Template{Idx{j}}]);
    end
else
    for j = 1:numel(Idx)
        T(j) = sum(Template(Idx{j}));
    end
end