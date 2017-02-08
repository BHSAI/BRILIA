function [Gx,Ax] = plotVvsCDR3mut(Output,varargin)
P = inputParser;
addRequired(P,'Output',@isnumeric);
addOptional(P,'CapValue',[],@isnumeric);
parse(P,Output,varargin{:});

Output = P.Results.Output;
CapValue = P.Results.CapValue;

%Select the inputs
Vmut = Output(:,1);
CDR3mut = Output(:,2);

%Generating the 2D histogram
Incr = 0.015;
Xval = 0:Incr:0.15;
Yval = 0:Incr:0.15;
PDFmat = zeros(length(Yval),length(Xval));
for r = 1:length(Yval)-1
    for c = 1:length(Xval)-1  
        PDFmat(r,c) = sum(Vmut>=Xval(c) & Vmut<Xval(c+1) & CDR3mut>=Yval(r) & CDR3mut<Yval(r+1));
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

%Formatting axes labels
set(Ax,'FontName','Arial','FontSize',12);
set(Ax,'Ydir','normal');
set(Ax,'box','on');
set(get(Ax,'XLabel'),'String','Vframe Mut Fraction','FontName','Arial','FontSize',14)
set(get(Ax,'YLabel'),'String','CDR3 Mut Fraction','FontName','Arial','FontSize',14)

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

[~, Ax] = formatPlotSize(Ax,3.3,3.3);
