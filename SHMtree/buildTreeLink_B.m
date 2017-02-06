%buildTreeLink_B perform the tree linking given a defined cluster, where
%the first seq is the root seq. 

%The way it works is this: 
%1) First seq "was" set as root, based on Template Count and Alignment Score
%2) Find the "next" seq with lowest score, and link it to root.
%3) Find the "next" seq with lowest score, and link it to root or other seq.
%4) repeat 3 until no more linking is left.
%* This linking step is indepedent of cutoff distance, as it will always
%link it.

function AncMap = buildTreeLink_B(Tdata,VDJheader)
H = getHeaderVar(VDJheader);

AncMap = [[1:size(Tdata,1)]' zeros(size(Tdata,1),4)]; %[ChildNum ParNum SHMHAMdist Template ClusterNum]
if H.TemplateLoc > 0
    AncMap(:,4) = cell2mat(Tdata(:,H.TemplateLoc));
    AncMap(isnan(AncMap(:,4)),4) = 1;
    AncMap(AncMap(:,4)==0,4) = 1;
else
    AncMap(:,4) = 1;
end

%Fill in the first entry
AncMap(1,3) = calcSHMHAMdist(Tdata{1,H.RefSeqLoc},Tdata{1,H.SeqLoc});

%Now link, starting from root and working down.
PairDist = calcPairDist(Tdata(:,H.SeqLoc),'shmham');
Active = ones(1,size(PairDist,2)) > 0;
Active(1) = 0;
while max(Active) == 1
    %Figure out active/inactive
    ActLoc = find(Active == 1);
    InactLoc = find(Active == 0);
    
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
