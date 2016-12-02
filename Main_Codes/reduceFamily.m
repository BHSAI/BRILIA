%reduceFamily(Xmatch) takes the data from findXmatch and filters out
%redundant gene returns, where X is either V, D, or J. 
%
%  VmatchR = reduceFamily(Vmatch) takes the Vmatch cell from findVmatch,
%  and then looks for matching family that have same LMRs. This redundancy
%  happens when the sample seq is too small, hence, there can be multiple
%  same matches. VmatchR will auto align all matches, find consensus, then
%  reduce the Vmatch into only unique families.

function UnqXmatch = reduceFamily(Xmatch)
%Sort rows based on gene name
Xmatch = sortrows(Xmatch,2);

%Determine all family names
Family = cell(size(Xmatch,1),1);
for j = 1:size(Xmatch,1)
    Name = Xmatch{j,2};
    DashLoc = strfind(Name,'-');
    if isempty(DashLoc)
        DashLoc = strfind(Name,'*');
    end
    if isempty(DashLoc)
        DashLoc = length(Name);
    end
    Family{j} = Name(1:DashLoc-1);
end
[~, ~, Map1] = unique(Family);

%Find the unique LMRs associated with each entry
LMR = cell2mat(Xmatch(:,3));
[~, ~, Map2] = unique(LMR,'rows');

%Now find the unique Family AND LMR matches
Map = [Map1 Map2];
[UnqMap, ~, ComboLoc] = unique(Map,'rows');

%Reduce the Xmatch cell based on family and LMR similarities.
UnqXmatch = cell(size(UnqMap,1),size(Xmatch,2));
for j = 1:max(ComboLoc)
    j1 = find(ComboLoc == j);
    if length(j1) > 1; %Condense similar entries
        Families = [Family{j1(1)} ': ' Xmatch{j1(1),2}];
        FamilyNum = [Xmatch{j1(1),1}];
        for k = 2:length(j1) %Concantenate Family name and # in file.
           Families = [Families '|' Xmatch{j1(k),2}];
           FamilyNum = [FamilyNum Xmatch{j1(k),1}];
        end
        UnqXmatch(j,3:end) = deal(Xmatch(j1(1),3:end));    
        UnqXmatch{j,2} = Families;
        UnqXmatch{j,1} = FamilyNum;
    else
        UnqXmatch(j,:) = deal(Xmatch(j1,:));    
    end
end