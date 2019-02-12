%makeConfettiMap will make a confetti map given the sizes of each cluster.
%
%  Confetti = makeConfettiMap(Sizes)
%
%  INPUT
%    GrpSizes: N-elem array of N cluster sizes
%    
%  OUTPUT
%    Confetti: sqrt(N) x sqrt(N) size square matrix map of clusters

function Confetti = makeConfettiMap(GrpSizes)

GrpSizes = GrpSizes(:);
BoxArea  = sum(GrpSizes);
BoxSide  = ceil(sqrt(BoxArea));
NormGrpSizes = sort(ceil(BoxSide^2/BoxArea * GrpSizes), 'descend');
NormBoxSide  = ceil(sqrt(sum(NormGrpSizes)));

Confetti = zeros(NormBoxSide, 'uint32');
for Num = 1:length(NormGrpSizes)
%     if mod(Num, round(length(NormGrpSizes)/10)) == 0
%         disp(Num/length(NormGrpSizes));
%     end
    [Rc, Cc] = find(Confetti == 0, 1);
    if isempty(Rc); break; end
    [~, ~, ~, Confetti] = findGroupCoord(Confetti, Num, NormGrpSizes(Num), Rc, Cc);
end

%findGroupCoord will find the map coordinates for the confetti map per
%region
function [R, C, Idx, GridOut] = findGroupCoord(Grid, Num, Size, Rc, Cc)
if Size <= 0 || Rc <= 0 || Cc <= 0 || Grid(Rc, Cc) ~= 0
    R = 0;
    C = 0;
    Idx = 0;
    return
end

R = zeros(floor(Size), 1);
C = zeros(floor(Size), 1);

%The next block is always the one that'll give you the most connectivity
%with existing block, going either clockwise or counterclockwise

%Initialize
s = 1;
R(s) = Rc;
C(s) = Cc;
Grid(Rc, Cc) = Num;

Ra = Rc; %active point
Ca = Cc; %active point
while s <= Size
    [Rn, Cn, IdxN] = findPerimeterMEX(size(Grid), Ra, Ca);
    Valids = Grid(IdxN) == 0;
    Rn = Rn(Valids);
    Cn = Cn(Valids);
    
    %Possible that you reached a dead end in space. Go back to another area
    if isempty(Rn)
        for k = 1:s-1
            [Rn, Cn, IdxN] = findPerimeterMEX(size(Grid), R(k), C(k)); 
            Valids = Grid(IdxN) == 0;
            Rn = Rn(Valids);
            Cn = Cn(Valids);
            if ~isempty(Rn); break; end
        end
    end
    
    %If it's still empty, find any empty spot closest to the cluster
    if isempty(Rn)
        [Rn, Cn] = find(Grid == 0);
    end
    if isempty(Rn); break; end
    
    AdjCount = zeros(length(Rn), 1);
    for k = 1:length(Rn)
        [~, ~, IdxO] = findPerimeterMEX(size(Grid), Rn(k), Cn(k)); 
        AdjCount(k) = sum(Grid(IdxO) == Num);
    end
    BestLoc = find(AdjCount == max(AdjCount));
    if length(BestLoc) > 1 %Find the closer one to Rc, Cc
        Dist = (Rc - Rn(BestLoc)).^2 + (Cc- Cn(BestLoc)).^2;
        BestDist = find(Dist == min(Dist));
        BestLoc = BestLoc(BestDist(1));
    end
    
    BestLoc = BestLoc(1);
    R(s) = Rn(BestLoc);
    C(s) = Cn(BestLoc);
    Grid(R(s), C(s)) = Num;
    Ra = R(s);
    Ca = C(s);
    s = s+1;
end

if nargout >= 3
    ZeroLoc = R == 0 | C == 0;
    Idx = sub2ind(size(Grid), R(~ZeroLoc), C(~ZeroLoc));
    if nargout >= 4
        GridOut = Grid;
    end
end