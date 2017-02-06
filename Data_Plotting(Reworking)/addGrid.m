%will add Grids to the images figures made in Matlab.


function addGrid(gca1)
hold(gca1,'on')
Xticks = get(gca1,'XTick');
dX = Xticks(2)-Xticks(1);
Xgrid = [Xticks(1)-dX/2:dX:Xticks(end)+dX/2];
Yticks = get(gca1,'YTick');
dY = Yticks(2)-Yticks(1);
Ygrid = [Yticks(1)-dY/2:dY:Yticks(end)+dY/2];

%Draw the Y grid lines
for i = 1:length(Xgrid)
   plot(gca1,[Xgrid(i) Xgrid(i)],[Ygrid(1) Ygrid(end)],'k')
   hold on    
end

%Draw the X grid lines
for i = 1:length(Ygrid)
   plot(gca1,[Xgrid(1) Xgrid(end)],[Ygrid(i) Ygrid(i)],'k')
   hold on    
end

hold off
