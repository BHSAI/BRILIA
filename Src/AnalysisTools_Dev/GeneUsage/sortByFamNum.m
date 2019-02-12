%sortByFamNum will take a list of gene names are reorganize by the gene
%family number, returning the resorted list and the index.

%Obtains the index to sort the List by gene family number
function [List, SortIdx] = sortByFamNum(List)
SortCell = cell(length(List), 1);
for k = 1:length(List)
    Name = parseGeneName(List{k});
    Num = convStr2Num(Name{1}(5:end));
    SortCell{k, 1} = [Name{1}(3) sprintf('%04.0f', Num)]; %The leading 0's are to help ensure correct sort
end
[~, SortIdx] = sort(SortCell);
List = List(SortIdx);