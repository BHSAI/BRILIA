function Gx = plotMutationCorr(varargin)
%See if user gave VDJdata and NewHeader, or MutData, or nothing
if isempty(varargin) || (~isempty(varargin) && isempty(varargin{1})) %Need to find file
    MutData = getMutationData;
elseif ~isempty(varargin)
    if ischar(varargin{1}) %Filename was given
        MutData = getMutationData(varargin{1});
        varargin(1) = [];
    elseif length(varargin) == 2 && iscell(varargin{1}) && iscell(varargin{2}) %VDJdata and VDJheader was given
        MutData = getMutationData(varargin{1:2});
        varargin(1:2) = [];
    end
else
    error('getMutationData: Check the inputs');
end
H = getHeaderVar(VDJheader);

%What remain of varargin, extract setting
P = inputParser;
addOptional(P,'CapValue',[],@isnumeric);
parse(P,varargin{:});
CapValue = P.Results.CapValue;

%Calculating the mutation frequencies
DelThis = sum(Output(:,1:4),2) == 0; %Just in case there are no mut in seq
Output(DelThis,:) = [];
Amut = Output(:,1); %Y axis
Cmut = Output(:,2); %X axis
Gmut = Output(:,3); %X axis
Tmut = Output(:,4); %Y axis

Atot = Output(:,5); %Y axis
Ctot = Output(:,6); %X axis
Gtot = Output(:,7); %X axis
Ttot = Output(:,8); %Y axis

ATmut = (Amut+Tmut)./(Atot+Ttot); %Normalized AT mutation
CGmut = (Cmut+Gmut)./(Ctot+Gtot); %Normalized CG mutation

%Generating the 2D PDF/histogram
Incr = 0.02;
Xval = 0:Incr:0.2;
Yval = 0:Incr:0.2;
PDFmat = zeros(length(Xval));
for r = 1:length(Yval)-1
    for c = 1:length(Xval)-1  
        PDFmat(r,c) = sum(CGmut>=Xval(c) & CGmut<Xval(c+1) & ATmut>=Yval(r) & ATmut<Yval(r+1));
    end
end

%Want to ensure that you can see most of the plots. So set a cap.
PDFvals = sort(PDFmat(:),'descend');
if isempty(CapValue)
    CapValue = round(PDFvals(2),-2);
end
PDFmat(PDFmat>CapValue) = CapValue;

%Setup the plot
[Gx, Ax] = formatPlotSize([],3.3,3.3);
Ix = image(PDFmat);
Ix.CDataMapping = 'scaled';
Cgray = colormap('gray');
Cgray = flipud(Cgray); %inverting color sheme;
colormap(Ax,Cgray);

%Draw the colorbar
c1 = colorbar;
c1.FontSize = 12;   
c1.Location = 'east';
ColorBarLim = get(c1,'Limits');
Round10 = ceil(CapValue^(1/10));
ColorBarTicks = [0:10^Round10:ColorBarLim(2)];

%Set the colorbar label to frequency count
ColorBarSpacing = ColorBarTicks/ColorBarLim(2)*max(PDFmat(:));
ColorBarTickLabel = cell(size(ColorBarSpacing));
for j = 1:length(ColorBarSpacing)
    if j < length(ColorBarSpacing)
        ColorBarTickLabel{j} = sprintf('%0.0f',ColorBarSpacing(j));
    else
        ColorBarTickLabel{j} = sprintf('>%0.0f',ColorBarSpacing(j));
    end
end
set(c1,'Ticks',ColorBarTicks,'TickLabels',ColorBarTickLabel);% ColorBarTicks = [0:0.1:1]*ColorBarLim(2);

%Fix the axis labels
% Xx = xlabel('(F_{C\rightarrowT} + F_{G\rightarrowA}) / (C + G)');
% Yx = ylabel('(F_{A\rightarrowG} + F_{T\rightarrowC}) / (A + T)');
Xx = xlabel('CG_{mut} / CG_{tot}');
Yx = ylabel('AT_{mut} / AT_{tot}');
set(Ax,'FontName','Arial','FontSize',14);
set(Xx,'FontName','Arial','FontSize',14);
set(Yx,'FontName','Arial','FontSize',14);
set(Ax,'Ydir','normal');

%Set the X and Y limits to match with the bin X and Y
Xlab = cell(1,length(Xval));
for k = 1:length(Xlab)
    Xlab{k} = sprintf('%0.2f',Xval(k));
end
Xtick = [0.5:length(Xval)+0.5];

Ylab = cell(1,length(Yval));
for k = 1:length(Ylab)
    Ylab{k} = sprintf('%0.2f',Yval(k));
end
Ytick = [0.5:length(Yval)+0.5];

%Taking every other label
Xtick = Xtick(1:2:end);
Ytick = Ytick(1:2:end);
Xlab = Xlab(1:2:end);
Xlab{1} = '0';
Ylab = Ylab(1:2:end);
Ylab{1} = '0';
set(Ax,'XTick',Xtick,'XTickLabel',Xlab,'YTIck',Ytick,'YTickLabel',Ylab);

[~, Ax] = formatPlotSize(Ax,2.7,2.7,'inner','AddBorder',0.05);

%==========================================================================
%Identify the 90% containment circle
Rad = [0:0.002:1];
Ct90 = length(ATmut)*0.90; %sum(PDFmat(:))*0.90;
AllDist = sqrt(ATmut.^2 + CGmut.^2);
for j = 1:length(Rad)
    Ct = sum(AllDist <= Rad(j));
    if Ct > Ct90;
        break
    else
        R90 = Rad(j);
    end
end
R90c = R90/Incr; %need to convert to axis values

%Draw the circle on a scatter plot
hold on
Xcir = [0:R90c/100:R90c];
Ycir = sqrt(R90c.^2-Xcir.^2);
plot(Ax,Xcir,Ycir,'--r','LineWidth',2);
%text(1,10,['Radius is ' num2str(R90,'%0.2f')],'FontSize',14,'Color',[1 0 0],'FontWeight','Bold');
hold off


% 
% 
% % 
% % ColorBarSpacing = ColorBarTicks;
% % ColorBarTickLabel = cell(size(ColorBarSpacing));
% % for j = 1:length(ColorBarSpacing)
% %     ColorBarTickLabel{j} = sprintf('%0.0f',ColorBarSpacing(j));
% % end
% % set(c1,'Ticks',ColorBarTicks,'TickLabels',ColorBarTickLabel);
% 
% 
% 
% 
% 
% 
% 
% for j = 1:length(ATmut)
%     Row = find(ATmut(j) < Xval);
%     if isempty(Row)
%         Row = 1;
%     else
%         Row = Row(1)-1;
%     end
%     Col = find(CGmut(j) < Yval);
%     if isempty(Col)
%         Col = 1;
%     else
%         Col = Col(1)-1;
%     end
%     PDFmat(Col,Row) = PDFmat(Col,Row) + 1; %Col and Row are flipped, since the "X" value must go by column, and the "Y" vlaue must go by row. 
% end
% 
% [~, Ax] = formatPlotSize([],3.3,3.3);
% scatter(Ax,CGmut,ATmut,'kx')



% %Chanage the tick range and labels
% XTicks = [1:2:length(Xval)];
% YTicks = [1:2:length(Yval)];
% set(Ax,'XLim',[1-0.5 length(Xval)+0.5],'XTick',XTicks)
% set(Ax,'YLim',[1-0.5 length(Yval)+0.5],'YTIck',YTicks)
% XTickLab = cell(1,length(XTicks));
% for k = 1:length(XTicks)
%     XTickLab{k} = sprintf('%0.2f',Xval(XTicks(k)));
% end
% YTickLab = cell(1,length(YTicks));
% for k = 1:length(YTicks)
%     YTickLab{k} = sprintf('%0.2f',Yval(YTicks(k)));
% end
% set(Ax,'XTickLabel',XTickLab,'YTickLabel',YTickLab)


