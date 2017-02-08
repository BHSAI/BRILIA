%plotGeneFamilyCombo will plot a VDJ gene usage scatter plot instead of
%histograms. 
function Gx = plotGeneFamilyCombo(VDJcombo,AxisLabels,varargin)
CustomColor = [1:-0.05:0]';
CustomColor = [CustomColor CustomColor ones(size(CustomColor,1),1)];

%DJ combination
DJcombo = squeeze(sum(VDJcombo,1));
UsedJ1 = sum(DJcombo,1) > 0;1
UsedD1 = sum(DJcombo,2) > 0;

%VD combination
VDcombo = squeeze(sum(VDJcombo,3));
UsedD2 = sum(VDcombo,1) > 0;
UsedV2 = sum(VDcombo,2) > 0;

%Do not remove unused stuff
UseAll = 1;
if UseAll == 1
    UsedJ1 = ones(size(UsedJ1))==1;
    UsedD1 = ones(size(UsedD1))==1;
    UsedV2 = ones(size(UsedV2))==1;
    UsedD2 = ones(size(UsedD2))==1;
end

DJcombo(:,UsedJ1 == 0) = [];
DJcombo(UsedD1 == 0,:) = [];
VDcombo(:,UsedD2 == 0) = [];
VDcombo(UsedV2 == 0,:) = [];

VD_DJcombo = [VDcombo' zeros(size(DJcombo,1),1) DJcombo];
Xlabels = [AxisLabels{1}(UsedV2); {''}; AxisLabels{3}(UsedJ1)];
Ylabels = [AxisLabels{2}(UsedD2)];

%Calculate the V and J total frequencies
VJtotFreq = sum(VD_DJcombo,1);
DtotFreq = sum(VD_DJcombo,2);

%Plot the VJ freq bars
Gx = figure;
set(Gx,'units','inches');
set(Gx,'position',[0 0 5 5]);
set(Gx,'PaperPosition',[0 0 6.5 size(VD_DJcombo,1)*6.5/(size(VD_DJcombo,2))]);




subplot(2,2,1)
bar(VJtotFreq)

%Plot the D freq bars
subplot(2,2,4)
barh(DtotFreq)

%Plot the scatter plot freq
subplot(2,2,3)
AX = gca;

%Plot the data
VDJidx = find(VD_DJcombo>0);
[Rx, Cx] = ind2sub(size(VD_DJcombo),VDJidx);
Fx = VD_DJcombo(VDJidx);
ScaleFx = 500/max(Fx);

scatter(Cx,Rx,Fx*ScaleFx,'fill','MarkerFaceColor','k');
colormap(CustomColor)
set(AX,'XTick',[1:length(Xlabels)],'YTick',[1:length(Ylabels)]);
set(AX,'XTickLabel',Xlabels,'YTickLabel',Ylabels,'XTickLabelRotation',90);
set(AX,'FontName','Arial','FontSize',12)
set(AX,'TickLength',[0 0])

%Fix the grids and plot sizes for tight fit
% addGrid(AX)
xlim([0.5 size(VD_DJcombo,2)+0.51])
ylim([0.5 size(VD_DJcombo,1)+0.51])
set(AX,'box','on')

%Extend the plots to fill space
set(AX,'OuterPosition',[0 0 1 1]);
TightInsets = get(AX,'TightInset');
Positions = [TightInsets(1) TightInsets(2) 1-TightInsets(3)-TightInsets(1) 1-TightInsets(4)-TightInsets(2)];
set(AX,'Position',Positions)

% 
% 
% 
% %Plot the images
% Gx{1} = figure;
% %image(DJcombo,'CDataMapping','scaled');
% DJidx = find(DJcombo>0);
% [Rx, Cx] = ind2sub(size(DJcombo),DJidx);
% Fx = DJcombo(DJidx);
% Sx = scatter(Cx,Rx,Fx*100,'fill','MarkerFaceColor','k')
% colormap(CustomColor)
% Ax{1} = gca;
% set(gca,'XTick',[1:size(DJcombo,2)],'YTick',[1:size(DJcombo,1)]);
% set(gca,'XTickLabel',AxisLabels{3}(UsedJ1),'YTickLabel',AxisLabels{2}(UsedD1),'XTickLabelRotation',90);
% set(gca,'FontName','Arial','FontSize',12)
% set(gca,'OuterPosition',[0 0 1 1])
% set(gca,'TickLength',[0 0])
% addGrid(gca)
% 
% Gx{2} = figure;
% VDcombo = VDcombo'; %Going to make Y axis the D's
% image(VDcombo,'CDataMapping','scaled');
% Ax{2} = gca;
% colormap(CustomColor)
% set(gca,'XTick',[1:size(VDcombo,2)],'YTick',[1:size(VDcombo,1)]);
% set(gca,'XTickLabel',AxisLabels{1}(UsedV2),'YTickLabel',AxisLabels{2}(UsedD2),'XTickLabelRotation',90);
% set(gca,'FontName','Arial','FontSize',12)
% set(gca,'OuterPosition',[0 0 1 1])
% set(gca,'TickLength',[0 0])
% addGrid(gca)
% 
% Gx{3} = figure;
% image(VJcombo,'CDataMapping','scaled');
% Ax{3} = gca;
% colormap(CustomColor)
% set(gca,'XTick',[1:size(VJcombo,2)],'YTick',[1:size(VJcombo,1)]);
% set(gca,'XTickLabel',AxisLabels{3}(UsedJ3),'YTickLabel',AxisLabels{1}(UsedV3),'XTickLabelRotation',90);
% set(gca,'FontName','Arial','FontSize',12)
% set(gca,'OuterPosition',[0 0 1 1])
% set(gca,'TickLength',[0 0])
% addGrid(gca)
% 
% % %Scale the figure heights proportionally to the tallest vertical box
% % MaxHeight = 7; %inch
% % MaxCount  = max([size(VDcombo) size(VJcombo) size(DJcombo)]);
% % ScaleRatio = MaxHeight/MaxCount;
% % 
% % for q = 1:3
% %     NewHeight = length(get(Ax{q},'Ytick')) * ScaleRatio;
% %     set(Gx{q},'PaperSize',[7 NewHeight]);
% %     set(Gx{q},'PaperPosition',[0 0 7 NewHeight]);
% % end
% 
% %Scale the height of the plot to be same bar height, and scale the width to
% %be same width.
% 
% %Find largest Left and Bottom border, width and height of the plot areas
% LeftBorder = 0;
% BotBorder = 0;
% for q = 1:3
%     BoxPos = get(Ax{q},'TightInset');
%     if BoxPos(1) > LeftBorder
%         LeftBorder = BoxPos(1);
%     end
%     if BoxPos(2) > BotBorder
%         BotBorder = BoxPos(2);
%     end
% end
% BoxPosition = [LeftBorder BotBorder 0.99-LeftBorder 0.99-BotBorder];
% for q = 1:3 %Standadize all positions
%     set(Ax{q},'Position',BoxPosition);
% end
