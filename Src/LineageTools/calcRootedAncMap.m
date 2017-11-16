%calcRootedAncMap perform the tree linking given a defined cluster, where
%the first seq is the most ancestral seq.
%  
%The way it works is this: 
%1) First seq "was" set as ancestor before puting the input in.
%2) Find the "next" seq with lowest score, and link it to root.
%3) Find the "next" seq with lowest score, and link it to root or other seq.
%4) repeat 3 until no more linking is left.
%
%  AncMap = calcRootedAncMap(PairDist)
%
%  INPUT
%    PairDist: a MxM matrix of distances between sequences, where the row
%      is the parent and col is the child sequences. 1st Row and Col is for
%      the root sequence.
%
%  OUTPUT
%    AncMap: a Mx3 matrix storing information about each child object's
%      parent object and the SHM distance from parent to child object.
%
%  NOTE
%    This does not have a cutoff distance as it will link everything to the
%    ancestral sequence.
%    
%    Unlike calcAncMap, calcRootedAncMap does not generate cyclic
%    dependencies since the root is defined and this prevents the issue.
function AncMap = calcRootedAncMap(PairDist)
%Ensure PairDist is a square matrix
[M, N] = size(PairDist);
if M ~= N %Check for linear form
    M = ceil(sqrt(2*N));
    if M*(M-1)/2 == N
        PairDist = squareform(PairDist);
    else
        error('%s: Input not a square distance matrix.', mfilename);
    end
end
PairDist(eye(size(PairDist))>0) = Inf; %Set diagonal to max to prevent self matching

%Now link, starting from root and working down.
AncMap = [[1:size(PairDist, 1)]' zeros(size(PairDist, 1), 2)]; %[ChildNum ParNum SHMHAMdist]
Active = ones(1, size(PairDist, 2), 'logical');
Active(1) = 0;
while max(Active) == 1
    %Figure out active/inactive
    ActLoc = find(Active);
    InactLoc = find(~Active);
    
    %Extract just the distance for active child, inactive parent
    PairDistT = PairDist(InactLoc, :);
    PairDistT = PairDistT(:, ActLoc);
    
    %Identify the next immediate child for each parent
    [ParentIdx, ChildIdx] = find(PairDistT == min(PairDistT(:)));
    ParentLoc = InactLoc(ParentIdx);
    ChildLoc = ActLoc(ChildIdx);
    
    %Ties are automatically broken. Ex: if 1 child belongs to 2 parents, 
    %based on the for loop below, the last parent wins.
    
    %For each child, fill in the parent
    for c = 1:length(ChildLoc)
        AncMap(ChildLoc(c), 2) = ParentLoc(c);
        AncMap(ChildLoc(c), 3) = double(PairDist(ParentLoc(c), ChildLoc(c)));
        Active(ChildLoc(c)) = 0;
    end
end
