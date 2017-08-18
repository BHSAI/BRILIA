%DEPRACATED, replaced with calcRootedAncMap
%
%buildTreeLink_B perform the tree linking given a defined cluster, where
%the first seq is the root seq. 

%The way it works is this: 
%1) First seq "was" set as root, based on Template Count and Alignment Score
%2) Find the "next" seq with lowest score, and link it to root.
%3) Find the "next" seq with lowest score, and link it to root or other seq.
%4) repeat 3 until no more linking is left.
%* This linking step is indepedent of cutoff distance, as it will always
%link it.

function AncMap = buildTreeLink_B(Seq,RefSeq) %Tdata,VDJheader)
error('%s: This is depracated. Replaced with calcRootedAncMap',mfilename);

%Ensure Seq is a cell of string, RefSeq is just a characterter
if ischar(Seq)
    Seq = {Seq};
end
if iscell(RefSeq)
    RefSeq = RefSeq{1};
end

%Fill in the first entry
AncMap = [[1:size(Seq,1)]' zeros(size(Seq,1),2)]; %[ChildNum ParNum SHMHAMdist]
AncMap(1,3) = calcSHMHAMdist(Seq{1},RefSeq);

%Now link, starting from root and working down.
PairDist = calcPairDist(Seq,'shmham');
Active = ones(1,size(PairDist,2),'logical');
Active(1) = 0;
while max(Active) == 1
    %Figure out active/inactive
    ActLoc = find(Active);
    InactLoc = find(~Active);
    
    %Extract just the distance for active child, inactive parent
    PairDistT = PairDist(InactLoc,:);
    PairDistT = PairDistT(:,ActLoc);
    
    %Identify the next immediate child for each parent
    [ParentIdx,ChildIdx] = find(PairDistT == min(PairDistT(:)));
    ParentLoc = InactLoc(ParentIdx);
    ChildLoc = ActLoc(ChildIdx);
    
    %Ties are automatically broken. Ex: if 1 child belongs to 2 parents,
    %based on the for loop below, the last parent wins.
    
    %For each child, fill in the parent
    for c = 1:length(ChildLoc)
        AncMap(ChildLoc(c),2) = ParentLoc(c);
        AncMap(ChildLoc(c),3) = double(PairDist(ParentLoc(c),ChildLoc(c)))/2; %Remember, calcPairDist doubled the value 
        Active(ChildLoc(c)) = 0;
    end
end
