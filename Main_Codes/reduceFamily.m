%reduceFamily is used by findGeneMatch to reduce VDJ matches with the same
%gene family, gene deletion, and alignment score. This is used to "reduce"
%the number of degenerate match outputs by condensing the family names into
%something like: "IGHV01: IGHV01-01*01, IGHV01-01*02". Will return all
%unique gene families. Selecting which one of the degenerate solution is
%best must be determined by a later function or by findGeneMatch.
%
%  ReducedMatch = reduceFamily(GeneMatch)
%
%  See also findGeneMatch, findVDJmatch

function UnqGeneMatch = reduceFamily(GeneMatch)
%Sort rows based on gene name
GeneMatch = sortrows(GeneMatch,2);

%Determine all family names
Family = cell(size(GeneMatch,1),1);
for j = 1:size(GeneMatch,1)
    Name = GeneMatch{j,2};
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
LMR = cell2mat(GeneMatch(:,3));
[~, ~, Map2] = unique(LMR,'rows');

%Now find the unique Family AND LMR matches
Map = [Map1 Map2];
[UnqMap, ~, ComboLoc] = unique(Map,'rows');

%Reduce the Xmatch cell based on family and LMR similarities.
UnqGeneMatch = cell(size(UnqMap,1),size(GeneMatch,2));
for j = 1:max(ComboLoc)
    q = find(ComboLoc == j);

    %Condense similar entries
    if length(q) > 1;
        %Concantenate Family name and # in file.
        Families = [Family{q(1)} ': ' GeneMatch{q(1),2}];
        FamilyNum = [GeneMatch{q(1),1}];
        for k = 2:length(q) 
           Families = [Families '|' GeneMatch{q(k),2}];
           FamilyNum = [FamilyNum GeneMatch{q(k),1}];
        end
        
        %Save the changes
        UnqGeneMatch(j,3:end) = deal(GeneMatch(q(1),3:end));    
        UnqGeneMatch{j,2} = Families;
        UnqGeneMatch{j,1} = FamilyNum;
    
    %No degenerate match.
    else
        UnqGeneMatch(j,:) = deal(GeneMatch(q,:));    
    end
end
