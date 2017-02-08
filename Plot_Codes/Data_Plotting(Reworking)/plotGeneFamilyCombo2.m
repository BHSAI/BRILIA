%plotGeneFamilyCombo will plot a VDJ gene usage scatter plot instead of
%histograms. 
function Gx = plotGeneFamilyCombo(VDJcombo,AxisLabels,varargin)
CustomColor = [1:-0.05:0]';
CustomColor = [CustomColor CustomColor ones(size(CustomColor,1),1)];

P = inputParser;
addParameter(P,'MaxDotSize',500,@isnumeric);
parse(P,varargin{:});
MaxDotSize = P.Results.MaxDotSize;

%DJ combination
DJcombo = squeeze(sum(VDJcombo,1));
UsedJ1 = sum(DJcombo,1) > 0;
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
DtotFreq = sum(VD_DJcombo,2)/2; %need to divide by 2, since you count twice (VD + DJ)

%Setup the figure sizes
Gx = figure;
set(Gx,'units','pixels');
set(Gx,'position',[50 50 600 600]);
set(Gx,'PaperPosition',[0 0 6.5 6.5]);

%Plot the VJ freq bars
subplot(2,2,1)
Ax1 = gca;
Bx1 = bar(Ax1,VJtotFreq);
hold(Ax1,'on')
plot(Ax1,size(VDcombo,1)*[1 1]+1,[0 2*max(VJtotFreq(:))],'-k','LineWidth',1)
hold(Ax1,'off')

%Plot the D freq bars
subplot(2,2,4)
Ax2 = gca;
Bx2 = barh(DtotFreq);

%Plot the scatter plot freq
subplot(2,2,3)
Ax3 = gca;
VDJidx = find(VD_DJcombo>0);
[Rx, Cx] = ind2sub(size(VD_DJcombo),VDJidx);
Fx = VD_DJcombo(VDJidx);
ScaleFx = MaxDotSize/max(Fx);

scatter(Cx,Rx,Fx*ScaleFx,'fill','MarkerFaceColor','k');
colormap(CustomColor)
set(Ax3,'XTick',[1:length(Xlabels)],'YTick',[1:length(Ylabels)]);
set(Ax3,'XTickLabel',Xlabels,'YTickLabel',Ylabels,'XTickLabelRotation',90);
set(Ax3,'FontName','Arial','FontSize',12)
set(Ax3,'TickLength',[0 0],'box','on')

hold(Ax3,'on')
plot(Ax3,size(VDcombo,1)*[1 1]+1,[0 2*size(VDcombo,2)],'-k','LineWidth',1)
hold(Ax3,'off')
%==========================================================================
%Begin formatting the plots to get the right positions, etc.

%1) Set the x limits
Xlimits = [0 size(VD_DJcombo,2)+1];
Ylimits = [0 size(VD_DJcombo,1)+1];
set(Ax1,'xlim',Xlimits,'xticklabels','','yticklabels','');
set(Ax2,'ylim',Ylimits,'xticklabels','','yticklabels','');
set(Ax3,'xlim',Xlimits,'ylim',Ylimits);

%2) set the positions
DL = 0.01; %Spacer

OutPos3 = [0 0 0.8 0.8];
set(Ax3,'OuterPosition',OutPos3)
TightPos = get(Ax3,'TightInset');
NewAxPos = [TightPos(1) TightPos(2) 0.8-TightPos(1) 0.8-TightPos(2)] ; %OutPos3(3)-sum(TightPos([1,3]))-0.03 OutPos3(4)-sum(TightPos([2,4]))-0.03];
set(Ax3,'Position',NewAxPos);

OutPos1 = [0 0.8 0.8 0.2];
NewAxPos1 = [NewAxPos(1) NewAxPos(2)+NewAxPos(4)+DL NewAxPos(3) 0.90*(1-NewAxPos(4)-NewAxPos(2))];
set(Ax1,'OuterPosition',OutPos1)
set(Ax1,'Position',NewAxPos1)

OutPos2 = [0.8 0 0.2 0.8];
NewAxPos2 = [NewAxPos(1)+NewAxPos(3)+DL NewAxPos(2) 0.90*(1-NewAxPos(3)-NewAxPos(1)) NewAxPos(4)];
set(Ax2,'OuterPosition',OutPos2)
set(Ax2,'Position',NewAxPos2)

%Redo the freq limits, based on maximum limits of all data.
AllMax = max([VJtotFreq(:); DtotFreq(:)]);
AllMax = AllMax*1.05;
if AllMax <= 0; return; end %****************NOTHING TO PLOT!
set(Ax1,'ylim',[0 AllMax]);
set(Ax2,'xlim',[0 AllMax]);

%Set the X and Y tick labels for frequency, based on maximums
YTick = [0 max(VJtotFreq(:))];
YLabels = {'0' num2str(max(VJtotFreq(:)),'%0.0f')};
set(Ax1,'YTick',YTick,'YTickLabel',YLabels);
XTick = [0 max(DtotFreq(:))];
XLabels = {'0' num2str(max(DtotFreq(:)),'%0.0f')};
set(Ax2,'XTick',XTick,'XTickLabel',XLabels,'XTickLabelRotation',90);

%Adjust the bar colors
Bx1.FaceColor = [0 0 0];
Bx2.FaceColor = [0 0 0];
