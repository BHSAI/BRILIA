%plotDensity will plot a density plot using images. Default resolution =
%500x500 pixels.

% X = [1 1 1 1 1 1 1 2 2 2 2 2]
% Y = [1 1 1 1 2 2 1 2 3 4 2 3]
% Xlim = [0 500]
% Ylim = [0 500]
function [Gx,Ax] = plotDensity(X,Y,Xlim,Ylim,DotRad)
%Determine the image size
Xlen = diff(Xlim);
Ylen = diff(Ylim);
Xscale = Xlen/500;
Yscale = Ylen/500;

%Determine the dotsize
[Xt,Yt] = meshgrid(-DotRad:DotRad,-DotRad:DotRad);
R = Xt.^2+Yt.^2;
%Gaus = exp(-R.^2/DotRad/10);
CircMask = R <= DotRad.^2;
[Rloc,Cloc] = find(CircMask == 1);
Rloc = Rloc - ceil(size(Xt,1)/2);
Cloc = Cloc - ceil(size(Xt,2)/2);

Xpx = round(X/Xscale);
Ypx = round(Y/Yscale);

%Plot the density circles
DensFig = zeros(500);
for j = 1:length(Xpx)
    Ridx = Ypx(j) + Rloc;
    Cidx = Xpx(j) + Cloc;
    DelThese = (Ridx <= 0) | (Ridx > size(DensFig,1)) | (Cidx <= 0) | (Cidx > size(DensFig,2));
    Ridx(DelThese) = [];
    Cidx(DelThese) = [];
    Rgaus = Rloc(DelThese == 0) + ceil(size(Xt,1)/2);
    Cgaus = Cloc(DelThese == 0) + ceil(size(Yt,1)/2);
    RCidx = sub2ind(size(DensFig),Ridx,Cidx);
    %GausIdx = sub2ind(size(CircMask),Rgaus,Cgaus);
    DensFig(RCidx) = DensFig(RCidx) + 1;
end

Gx = figure;
Ix = imagesc(DensFig);
Ax = get(Ix,'parent');
colormap('jet');

%Format Stuff
set(Ax,'Ydir','normal')
XTickVal = [0:Xlim(2)/5:Xlim(2)];
YTickVal = [0:Ylim(2)/5:Ylim(2)];
XTickLab = cell(1,length(XTickVal));
YTickLab = cell(1,length(YTickVal));
for k = 1:length(XTickVal)
    XTickLab{k} = sprintf('%0.2f',XTickVal(k));
    YTickLab{k} = sprintf('%0.2f',YTickVal(k));
end
PlotXTick = round(XTickVal/Xscale);
PlotYTick = round(YTickVal/Yscale);
PlotXTick(1) = 1;
PlotYTick(1) = 1;
set(Ax,'XTick',PlotXTick,'YTick',PlotYTick,'XTickLabel',XTickLab,'YTickLabel',YTickLab);

%Set the colorbar info
c1 = colorbar;
c1.FontSize = 12;   
c1.Location = 'east';
ColorBarLim = get(c1,'Limits');
ColorBarTicks = [0:ColorBarLim(2)/10:ColorBarLim(2)];

ColorBarSpacing = ColorBarTicks/ColorBarLim(2)*max(DensFig(:));
ColorBarTickLabel = cell(size(ColorBarSpacing));
for j = 1:length(ColorBarSpacing)
    ColorBarTickLabel{j} = sprintf('%0.0f',ColorBarSpacing(j));
end
set(c1,'Ticks',ColorBarTicks,'TickLabels',ColorBarTickLabel);
set(c1,'Color',[1 1 1]);
