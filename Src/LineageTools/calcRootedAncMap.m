%calcRootedAncMap perform the tree linking given a defined cluster, where
%the first seq is the most ancestral seq. The steps are as folllows:
%1) The PairDist matrix already has the 1st row & col as the "ancestor"
%2) Find the next seq with lowest score, and link it to root
%3) Repeat linking until all sequence is linked
%
%  AncMap = calcRootedAncMap(PairDist)
%
%  INPUT
%    PairDist: a MxM matrix of distances between sequences, where the row
%      is the parent and col is the child sequences. The First Row and Col
%      MUST be the root sequence.
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
%    dependencies since the root is defined, which prevents the issue.
%
function AncMap = calcRootedAncMap(PairDist)
[M, N] = size(PairDist);
if M ~= N
    PairDist = squareform(PairDist, 'tomatrix');
end
PairDist(1:size(PairDist, 1)+1:end) = Inf; %Prevent linking to self

%Now link starting from root and working down
AncMap = [[1:size(PairDist, 1)]' zeros(size(PairDist, 1), 2)]; %#ok<NBRAK> %[Child Par Dist]
Active = ones(1, size(PairDist, 1), 'logical');
Active(1) = 0;
ActiveNum = sum(Active);
while ActiveNum > 0
    %Identify the next immediate child for each parent
    InactiveIdx = find(~Active);
    ActiveIdx   = find( Active);
    PairDistT = PairDist(InactiveIdx, ActiveIdx);
    [R, C] = find(PairDistT == min(PairDistT(:)));
    RIdx = InactiveIdx(R);
    CIdx =   ActiveIdx(C);
    
   %Assign a parent to each child. If 1 child has > 1 parent, 1st parent wins.
    for q = 1:length(CIdx)
        if Active(CIdx(q)) %Still missing parent
            AncMap(CIdx(q), 2) = RIdx(q);
            AncMap(CIdx(q), 3) = PairDist(RIdx(q), CIdx(q));
            Active(CIdx(q)) = 0;
            ActiveNum = ActiveNum - 1;
        end
    end
end
