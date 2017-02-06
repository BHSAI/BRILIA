%plotTDTmatrix will plot the frequence of NT in the N regions
function Gx = plotTDTmatrix(TDTmatrix, TDTcomp, varargin)
NTlabel = {'A' 'C' 'G' 'T'};
FigFontSize = 9;
FigH = 2; %inches
FigW = 1.65; %inches
AddBorder = 0.05;

%Normaling input for % based
TDTmatrix = round(TDTmatrix/sum(TDTmatrix(:))*100,1);
if sum(TDTmatrix(:)) ~=  100 %rounding error. punt it to the highest value.
    [MaxRow, MaxCol] = find(TDTmatrix == max(TDTmatrix(:)));
    MaxRow = MaxRow(1);
    MaxCol = MaxCol(1);
    TDTmatrix(MaxRow,MaxCol) = (100 - sum(TDTmatrix(:))) + TDTmatrix(MaxRow,MaxCol);
end

%Formatting TDTcomp for 2 decimal points that round to 1.
TDTcomp = round(100*TDTcomp/sum(TDTcomp(:)));
if sum(TDTcomp) ~= 100
    MaxCol = find(TDTcomp == max(TDTcomp(:)));
    MaxCol = MaxCol(1);
    TDTcomp(MaxCol) = (100 - sum(TDTcomp)) + TDTcomp(MaxCol);
end

%Plot the bulk composition bar graph
Gx = figure;
subplot(1,2,1)
B1 = bar(1:4,TDTcomp);
B1.FaceColor = [0 0 0];
Ax1 = gca;

%Establish the XY ticks
Xlim = [0.5 4.5];
Xtick = [1:4];
Xlab = NTlabel;
if max(TDTcomp) > 100
    Ymax = (ceil(max(TDTcomp)/10)+1)*10;
else
    Ymax = 100;
end
Ylim = [0 Ymax];
Ytick = [0:10:Ymax];
set(Ax1,'FontName','Arial','FontSize',FigFontSize,'Xlim',Xlim,'Xtick',Xtick,'XTickLabel',Xlab,'Ylim',Ylim,'Ytick',Ytick,'YTickLabel','','TickDir','out');
set(get(Ax1,'title'),'string','nt comp.','FontName','Arial','FontSize',FigFontSize+3);

%Add numbers to the bars
for q = 1:4
    text(q,TDTcomp(q),num2str(TDTcomp(q),'%2.0f'),'vert','bot','FontSize',FigFontSize,'HorizontalAlignment','center');
end

%Add the number of flips here
if ~isempty(varargin)
    Left2Right = varargin{1};
    Yranges = get(gca,'Ylim');
    Ymax = Yranges(2);
    Xranges = get(gca,'Xlim');
    Xmid = sum(Xranges)/2;
    TextAdd = sprintf('Coding Side: %2.0f %s %2.0f%%',Left2Right(1)*100,char(177),Left2Right(3)*100);
    text(Xmid,0.95*Ymax,TextAdd,'FontSize',FigFontSize,'FontName','Arial','HorizontalAlignment','center');
end

set(Ax1,'OuterPosition',[0 0 0.5 1])

%--------------------------------------------------------------------------
%Draw the pairwise mutation matrix now
subplot(1,2,2)

image(TDTmatrix,'CDataMapping','scaled')
Ax2 = gca;
colormap(gray)
set(gca,'XTick',1:length(NTlabel))
set(gca,'XTickLabel',NTlabel)
set(gca,'XAxisLocation','top')
set(gca,'YTick',1:length(NTlabel))
set(gca,'YTickLabel',NTlabel)
set(gca,'YTickLabel',NTlabel)
set(gca,'FontSize',FigFontSize)

Crange = [0 50];
set(gca,'Clim',Crange);
set(get(Ax2,'title'),'string','paired nt comp.','FontName','Arial','FontSize',FigFontSize+3);

%Write the text
for r = 1:4
    for c = 1:4
        h1 = text(c,r,num2str(TDTmatrix(r,c),'%2.1f'),'VerticalAlignment','middle','FontSize',FigFontSize,'HorizontalAlignment','center');
        if TDTmatrix(r,c) > diff(Crange)/2
            h1.Color = [0 0 0];
        else
            h1.Color = [1 1 1];
        end
    end
end

set(Ax2,'OuterPosition',[0.5 0 0.5 1]);

%--------------------------------------------------------------------------
%Set the Ax positions to match with each other.
set(Gx,'units','inch');
FigPos = [0 0 FigW*2 FigH];
set(Gx,'Position',FigPos,'PaperPosition',FigPos,'PaperSize',FigPos(3:4));

%Change and for mat axes properties, AFTER rescaling figure. to take adv of
%normalized units.
set(Ax1,'units','inch');
set(Ax2,'units','inch');
OutPos1 = get(Ax1,'outerposition');
Tight1 = get(Ax1,'tightinset') + AddBorder;
AxPos1 = OutPos1 + [Tight1(1) Tight1(2) -sum(Tight1([1 3])) -sum(Tight1([2 4]))];

OutPos2 = get(Ax2,'outerposition');
Tight2 = get(Ax2,'tightinset') + AddBorder;
Tight2(2) = Tight1(2); %Setting the bottom to be the same
AxPos2 = OutPos2 + [Tight2(1) Tight2(2) -sum(Tight2([1 3])) -sum(Tight2([2 4]))];

set(Ax1,'position',AxPos1)
set(Ax2,'position',AxPos2)
