%calcAncMap will take a PairDist matrix, and then link the numbers using
%nearest neighbor method, with priority on shorter distances.
%
%  AncMap = calcAncMap(PairDist)
%
%  AncMap = calcAncMap(PairDist, TriangleSide)
%
%  INPUT
%    PairDist: a MxM dissimilarity matrix, or matrix storing a distance
%      value from paired objects. Diagonal values will not be used for
%      linking. 
%    TriangleSide ['botleft' 'topright' 'both']: Specify which side of
%      PairDist symmetric matrix to use for finding parent-distance pairs.
%      Defaults to 'both' in case there is an assymmetric distance matrix.
%
%  OUTPUT
%    AncMap: a Mx3 matrix storing information about each child object's
%      parent object and the distance from parent to child object.
%
%  NOTE
%    Each column is the child, and the row is the parent. So PairDist(1,2)
%    = 0.2 would mean Distance from Object 1 parent to Object 2 child is
%    0.2
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
function AncMap = calcAncMap(PairDist, varargin)
%Ensure PairDist is a square matrix
[M, N] = size(PairDist);
if M ~= N
    try
        PairDist = squareform(PairDist);
    catch
        error('%s : Input is not a square distance matrix', mfilename);
    end
end
PairDist(eye(size(PairDist))>0) = Inf; %Set diagonal Inf to prevent self matching

%Determine which triangle side to use
TriangleSide = 'both';
if ~isempty(varargin)
    TriangleSide = lower(varargin{1});
end

%Create linkup using shortest distance, with priority on short seq
UnqDist = unique(PairDist(:)); %Unique and sorted PairDist values
AncMap = [[1:size(PairDist, 1)]' zeros(size(PairDist, 1), 2)]; %#ok<NBRAK> %[SeqNum AncNum Distance]
Active = ones(size(PairDist, 1), 1, 'logical');                %Keeps track of active col/row
k = 1;
while any(Active) && k < length(UnqDist)  %it's k < length(UnqDist), and not "<=" to stop BEFORE linking Inf distances
    %Remove certain (R == C) pairs to prevent cyclic dep.
    [R, C] = find(PairDist == UnqDist(k));
    switch TriangleSide
        case 'botleft'
            KeepLoc = C < R;
            R = R(KeepLoc);
            C = C(KeepLoc);
        case 'topright'
            KeepLoc = R < C;
            R = R(KeepLoc);
            C = C(KeepLoc);
    end
        
    %Nearest neighbor linking, starting with shortest distance first
    for q = 1:length(R)
        if Active(C(q)) %Still missing parent
            AncMap(C(q), 2) = R(q);
            AncMap(C(q), 3) = PairDist(R(q), C(q));
            Active(C(q)) = 0;
        end
    end
    k = k+1;
end