%calcAncMap will take a PairDist matrix, and then link the numbers using
%nearest neighbor method, with priority on shorter distances.
%
%  AncMap = calcAncMap(PairDist)

function AncMap = calcAncMap(PairDist)
%Create linkup using shortest distance, with priority on short seq
AncMap = [[1:size(PairDist,1)]' zeros(size(PairDist,1),2) ones(size(PairDist,1),1)]; %[SeqNum AncNum Distance Frequency]
Active = ones(size(PairDist,1),1) == 1; %Keeps track of active col/row
UnqDist = unique(PairDist(:));
UnqDist(end) = [];
k = 1;
while sum(Active) > 0 && k <= length(UnqDist)   
    %Perform nearest neighbor linking, starting with shortest distance
    %first
    [R, C] = find(PairDist == UnqDist(k));    
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
