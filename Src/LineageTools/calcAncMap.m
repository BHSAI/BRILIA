%calcAncMap will take a PairDist matrix, and then link the numbers using
%nearest neighbor method, with priority on shorter distances.
%
%  AncMap = calcAncMap(PairDist)
%
%  AncMap = calcAncMap(PairDist, TriangleSide)
%
%  INPUT
%    PairDist: a MxM dissimilarity matrix, or matrix storing a distance
%      value from paired objects. 
%    TriangleSide ['lower' 'upper' 'both']: Specify which side of
%      PairDist symmetric matrix to use for finding parent-distance pairs.
%      Defaults to 'both' in case there is an assymmetric distance matrix.
%
%  OUTPUT
%    AncMap: a Mx3 matrix storing information about each child object's
%      parent object and the distance from parent to child object.
%
%  NOTE
%    In PairDist, each column is the child, and the row is the parent. So
%    PairDist(1,2) = 0.2 would mean Distance from Seq 1 to Seq 2 is 0.2
%
%    Unlike calcRootedAncMap, calcAncMap can generate cyclic dependencies
%    if given an asymmetric PairDist matrix with TriangleSide = 'both'.
%
%  EXAMPLE
%    PairDist = [
%         0.0000    6.5000   14.0000   29.0000    5.5000   12.5000;  
%         8.0000    0.0000    7.5000   20.5000   13.5000   22.5000; 
%        16.0000    8.0000    0.0000    7.0000   21.5000   30.5000; 
%        30.5000   20.5000    6.5000    0.0000   38.0000   51.0000; 
%         7.5000   14.0000   21.5000   38.5000    0.0000    5.0000;
%        15.5000   24.0000   31.5000   52.5000    6.0000    0.0000]; 
%
%    AncMap = calcAncMap(PairDist)
%    AncMap =
%         1.0000    5.0000    7.5000
%         2.0000    1.0000    6.5000
%         3.0000    4.0000    6.5000
%         4.0000    3.0000    7.0000
%         5.0000    1.0000    5.5000
%         6.0000    5.0000    5.0000
%
function AncMap = calcAncMap(PairDist, TriangleSide)
[M, N] = size(PairDist);
if M ~= N
    PairDist = squareform(PairDist, 'tomatrix');
end
PairDist(1:size(PairDist, 1)+1:end) = Inf; %Prevent linking to self

%Determine which triangle side to use
if nargin == 2
    if any(startsWith({'lower', 'botleft'}, TriangleSide, 'ignorecase', true))
        PairDist(triu(ones(size(PairDist), 'logical'),  1)) = Inf;
    elseif any(startsWith({'upper', 'topright'}, TriangleSide, 'ignorecase', true))
        PairDist(tril(ones(size(PairDist), 'logical'), -1)) = Inf;
    end
end

%Nearest neighbor linking, starting with shortest distance first
UnqDist = unique(PairDist(:));
AncMap = [[1:size(PairDist, 1)]' zeros(size(PairDist, 1), 2)]; %#ok<NBRAK> %[SeqNum AncNum Distance]
Active = ones(size(PairDist, 1), 1, 'logical');                %Keeps track of active col/row
ActiveNum = size(PairDist, 1);                                 %Keeps track of how many active col/row there are
k = 1;
while ActiveNum > 0 && k < length(UnqDist)  %it's k < length(UnqDist), and not "<=", to stop BEFORE linking Inf distances
    [RIdx, CIdx] = find(PairDist == UnqDist(k));
    for q = 1:length(CIdx)
        if Active(CIdx(q)) %Still missing parent
            AncMap(CIdx(q), 2) = RIdx(q);
            AncMap(CIdx(q), 3) = PairDist(RIdx(q), CIdx(q));
            Active(CIdx(q)) = 0;
            ActiveNum = ActiveNum - 1;
        end
    end
    k = k+1;
end