%calcAncMap will take a PairDist matrix, and then link the numbers using
%nearest neighbor method, with priority on shorter distances.
%
%  AncMap = calcAncMap(PairDist)
%
%  AncMap = calcAncMap(PairDist,TriangleSide)
%
%  INPUT
%    PairDist: a MxM dissimilarity matrix, or matrix storing a distance
%      value from paired objects. Diagonal values will not be used for
%      linking. 
%    TriangleSide ['botleft' 'topright']: Specify which side of PairDist
%      symmetric matrix to use for finding parent-distance pairs. Defaults
%      to 'both' in case there is an assymmetric distance matrix.
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
%    due if give Both triangle sides and an assymmetric matrix. 

function AncMap = calcAncMap(PairDist,varargin)
%Ensure PairDist is a square matrix
[M, N] = size(PairDist);
if M ~= N %Check for linear form
    M = ceil(sqrt(2*N));
    if M*(M-1)/2 == N
        PairDist = squareform(PairDist);
    else
        error('%s : Input not a square distance matrix');
    end
end
PairDist(eye(size(PairDist))>0) = Inf; %Set diagonal to max to prevent self matching

%Determine which triangle side to use
TriangleSide = 'both';
if ~isempty(varargin)
    TriangleSide = lower(varargin{1});
end

%Determine all unique distances
UnqDist = unique(PairDist(:)); %Unique and sorted PairDist values
UnqDist(end) = []; %Removing Inf as a unique PairDist.

%Create linkup using shortest distance, with priority on short seq
AncMap = [[1:size(PairDist,1)]' zeros(size(PairDist,1),1) zeros(size(PairDist,1),1)]; %[SeqNum AncNum Distance]
Active = ones(size(PairDist,1),1) == 1; %Keeps track of active col/row
k = 1;
while sum(Active) > 0 && k <= length(UnqDist)   
    %Find the shortest distance pairs
    [R, C] = find(PairDist == UnqDist(k));
    
    %Remove certain RC pairs on PairDist matrix to prevent cyclic dep.
    switch TriangleSide
        case 'botleft'
            DelIdx = C > R;
            R(DelIdx) = [];
            C(DelIdx) = [];
        case 'topright'
            DelIdx = R > C;
            R(DelIdx) = [];
            C(DelIdx) = [];
    end
        
    %Nearest neighbor linking, starting with shortest distance first
    for q = 1:length(R)
        Rc = R(q);
        Cc = C(q);
        if Active(Cc) == 1 %Still missing parent
            AncMap(Cc,2) = Rc;
            AncMap(Cc,3) = PairDist(Rc,Cc);
            Active(Cc) = 0;
        end
    end
    k = k+1;
end
