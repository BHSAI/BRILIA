%getGridCoord will find the Row and Column coordinates of a subplot in a
%figure, given all the positions of all the unique axes.
%
%  GridCoord = getGridCoord(AllPositions)
%
%  INPUT
%    AllPositions: Nx4 outer positions of N unique axes in a multi-plot
%    figure. Use getUniqueAxes.m . 
%
%  OUTPUT
%    GridCoord: Nx2 row vs column coordinate for each axes.

function GridCoord = getGridCoord(AllPositions)
GridCoord = zeros(size(AllPositions, 1), 2);
XPos = [AllPositions(:, 1) sum(AllPositions(:, [1 3]), 2)];
Active = ones(size(XPos, 1), 1, 'logical');
MinPos = min(XPos(:,2));
CurCol = 1;
while any(Active)
    CurRowIdx = XPos(:,1) <= MinPos & Active;
    GridCoord(CurRowIdx, 2) = CurCol;
    CurCol = CurCol + 1;
    Active(CurRowIdx) = 0;
    MinPos = min(XPos(Active, 2));
end

%Find the rows next
YPos = [AllPositions(:, 2) sum(AllPositions(:, [2 4]), 2)];
Active = ones(size(YPos, 1), 1, 'logical');
MaxPos = max(YPos(:, 1));
CurCol = 1;
while any(Active)
    CurColIdx = YPos(:, 2) >= MaxPos & Active;
    GridCoord(CurColIdx, 1) = CurCol;
    CurCol = CurCol + 1;
    Active(CurColIdx) = 0;
    MaxPos = max(YPos(Active, 1));
end

%Invert GridCoord so that row 1-N = bot-top, and col 1-N = left to right.
GridCoord(:, 1) = max(GridCoord(:, 1)) - GridCoord(:, 1) + 1;