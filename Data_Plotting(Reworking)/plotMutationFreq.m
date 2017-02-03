function Gx = plotMutationFreq(varargin)
NTlabel1 = {'A_0' 'C_0' 'G_0' 'T_0'}; %Just keep this as default.
NTlabel2 = {'A_1' 'C_1' 'G_1' 'T_1'}; %Just keep this as default.

%Look the "Legend" field and the data field
InputLoc = 1:length(varargin);
Legend = '';
NumPlotCol = 2; %Number of subplot columns
for k = 1:length(InputLoc)
    if ischar(varargin{k})
        if strcmpi(varargin{k},'Legend') %Title of each plot
            Legend = varargin{k+1};
            InputLoc(k:k+1) = [];
        end        
        if strcmpi(varargin{k},'PlotsPerRow') %Number of plots per row
            NumPlotCol = varargin{k+1};
        end
    end    
end

%Look for the data
AllMutMat = varargin(InputLoc);
if length(AllMutMat) < NumPlotCol
    NumPlotCol = length(AllMutMat);
end
NumPlotRow = ceil(length(AllMutMat)/NumPlotCol);

%Finalize Legend
if isempty(Legend)
    Legend = repmat({''},1,length(InputLoc));
end
if length(Legend) ~= length(AllMutMat)
    error('Error: Legend number not same as the input number')
end

%Plot each data
Gx = figure;
for q = 1:length(AllMutMat)
    %Extract the current mutation matrix and name
    XmatOrig = AllMutMat{q};
    Xmat = XmatOrig/sum(XmatOrig(:));

    if isempty(Legend);
        Xname = '';
    else
        Xname = Legend{q};
    end
    
    %Determine min max ranges
    Crange = [0 ceil(max(Xmat(:))/0.1)*0.1];
    if Crange(2) == 0 %Empty data
        Crange(2) = 1;
    end
    
    %Calculate the position based on the location and number of plots
    Pwidth = 1/NumPlotCol;
    Pheight = 1/NumPlotRow;
    [Cnum, Rnum] = ind2sub([NumPlotCol NumPlotRow],q); %This is inverted referencing because suplot goes left to right, whereas matrix are read top to bot.    
    Yloc = (NumPlotRow - Rnum)*Pheight;
    Xloc = (Cnum-1)*Pwidth;
    
    subplot(NumPlotRow,NumPlotCol,q)
    image(Xmat,'CDataMapping','scaled')
    set(gca,'XTick',1:4,'XTickLabel',NTlabel1,'XAxisLocation','top');
    set(gca,'YTick',1:4,'YTickLabel',NTlabel2);
    set(gca,'FontSize',24,'OuterPosition',[Xloc Yloc Pwidth Pheight]);
    set(get(gca,'Title'),'String',Xname,'FontSize',24);    

    %Setting the colors    
    colormap('jet')
    set(gca,'Clim',Crange);
    c1 = colorbar;
    c1.FontSize = 16;   
    
    ColorBarLim = get(c1,'Limits');
    ColorBarTicks = [0:0.05:ColorBarLim(2)];
    ColorBarTickLabel = cell(size(ColorBarTicks));
    for j = 1:length(ColorBarTicks)
        ColorBarTickLabel{j} = sprintf('%0.2f',ColorBarTicks(j));
    end
    set(c1,'Ticks',ColorBarTicks,'TickLabels',ColorBarTickLabel);
    
    %Write the text
    for r = 1:4
        for c = 1:4
            h1 = text(c,r,num2str(XmatOrig(r,c),'%0.0f'),'FontSize',16);
            h1.VerticalAlignment = 'middle';
            h1.HorizontalAlignment = 'center';
            if (Xmat(r,c) > diff(Crange)/4) && (Xmat(r,c) < 3/4*diff(Crange))
                h1.Color = [0 0 0];
            else
                h1.Color = [1 1 1];
            end
        end
    end
end