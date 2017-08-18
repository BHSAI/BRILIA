%function condenseGeneMap will convert a large DB into a smaller one for
%improved alignment speed. For instance, there is no need to do 400 V gene
%alignements if 200 of them are very similar to another one or are alleles.
%This cuts down the alignment time, though the actual V must be determine
%using the full sequence.

function Ymap = condenseGeneMap(Xmap, MapHeader)
disp('Under construction.')
return;
M = getMapHeaderVar(MapHeader);

%Get the gene alleles grouped up
GeneNames = cell(size(Xmap, M.GeneLoc), 3);
for j = 1:size(Xmap, 1)
    [Fam, Gene, Allele] = parseGeneName(Xmap{j, M.GeneLoc});
    GeneNames(j, :) = [Fam, Gene, Allele];
end

%Sort by unique FamNum and GeneNum
SortCell = cell(size(Xmap, 1), 1);
for j = 1:size(Xmap, 1)
    SortCell{j} = [GeneNames{j, 1} GeneNames{j, 2}];
end
[UnqSortCell, ~, UnqIdx] = unique(SortCell);

%For every unique FamNum - GeneNum, determine if alleles are similar enough
UnqIdx2 = zeros(size(UnqIdx));
for y = 1:max(UnqIdx)
    Idx = find(UnqIdx == y);
    if length(Idx) == 1; continue; end
    pause
    Seq = Xmap(Idx, M.SeqLoc)
    PDist = calcPairDist(Seq, 'identity', 2);
end


%Trim out the anchor dist to include 104C but nothing afterwards
DelLoc = zeros(size(Xmap, 1), 1, 'logical');
for j = 1:size(Xmap, 1)
    Xmap{j, M.MiscLoc} = j; %Stores the absolute index number
    AnchorPos = Xmap{j, M.AnchorLoc};
    if isempty(Xmap{j, M.SeqLoc})
        DelLoc(j) = 1;
        continue
    end
    if AnchorPos == 0
        continue
    else
        Xmap{j, M.SeqLoc} = Xmap{j, M.SeqLoc}(1:end-AnchorPos+3);
        Xmap{j, M.AnchorLoc} = 3;
    end
end

%Delete the empty ones
Xmap(DelLoc, :) = [];



