function varargout = plotGeneUsage(GeneUsage, varargin)
P = inputParser;
addParameter(P, 'FigWidth', 6, @(x) isnumeric(x) && x > 1);
addParameter(P, 'FigHeight', 6, @(x) isnumeric(x) && x > 1);
addParameter(P, 'YLim', [0 1.2], @(x) isnumeric(x) && length(x) == 2);
addParameter(P, 'XLim', [0.5 4.5], @(x) isnumeric(x) && length(x) == 2);
addParameter(P, 'XTickLabel', {'A', 'C', 'G', 'T'}, @(x) iscell(x) && length(x) == 4);
addParameter(P, 'FaceColor', [0 0 0], @(x) isnumeric(x) && length(x) == 3);
addParameter(P, 'FontName', 'Arial', @(x) ischar(x));
addParameter(P, 'FontSize', 12, @(x) isnumeric(x) && x > 1);
addParameter(P, 'SaveSubDir', 'Analysis', @ischar);
addParameter(P, 'DefaultSaveAs', 'HotMotifBarGraph.csv', @ishcar);
addParameter(P, 'Visible', 'on', @(x) ischar(x) && ismember(lower(x), {'on', 'off'}));
addParameter(P,'MaxDotSize', 500, @isnumeric);
[Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInput(P, varargin{:});
if ReturnThis
   varargout = {Ps, Pu, ExpPs, ExpPu};
   return;
end
FigHeight = Ps.FigHeight;
FigWidth = Ps.FigWidth;
XLim = Ps.XLim;
YLim = Ps.YLim;
XTickLabel = Ps.XTickLabel;
FaceColor = Ps.FaceColor;
FontName = Ps.FontName;
FontSize = Ps.FontSize;
SaveSubDir = Ps.SaveSubDir;
DefaultSaveAs = Ps.DefaultSaveAs;
Visible = Ps.Visible;
MaxDotSize = P.Results.MaxDotSize;

switch length(GeneUsage.Names)
    case 3
        Chain = 'H';
    case 2
        Chain = 'L';
    case 5
        Chain = 'HL';
    otherwise
        error('%s: Could not determine chain type.');
end

CustomColor = [1:-0.05:0]';
CustomColor = [CustomColor CustomColor ones(size(CustomColor,1),1)];

if contains(Chain, 'H')
	VDJcombo = GeneUsage.Count(:, :, 1:3);
    Vcount = sum(sum(VDJcombo, 2), 3);
    Jcount = sum(sum(VDJcombo, 1), 2);
    VDcount = sum(VDJcombo, 3)';
    DJcount = squeeze(sum(VDJcombo, 1));

    VDDJtotal = [VDcount zeros(size(DJcount,1),1) DJcount];
    VJtotal = [Vcount(:); 0; Jcount(:)];
    Dtotal = sum(sum(VDJcombo, 1), 3);
    Xlabels = [GeneUsage.Names{1}, {''}, GeneUsage.Names{3}];
    Ylabels = GeneUsage.Names{2};
end

Gx = figure('Visible', Visible);
set(Gx, 'units', 'pixels');
set(Gx, 'position', [50 50 600 600]);
set(Gx, 'PaperPosition', [0 0 6.5 6.5]);

%Plot the VJ freq bars
subplot(2,2,1)
Ax1 = gca;
Bx1 = bar(Ax1, VJtotal);
hold(Ax1,'on')
plot(Ax1,size(VDcount,2)*[1 1]+1,[0 2*max(VJtotal(:))],'-k','LineWidth',1)
hold(Ax1,'off')

%Plot the D freq bars
subplot(2,2,4)
Ax2 = gca;
Bx2 = barh(Dtotal);

%Plot the scatter plot freq
subplot(2,2,3)
Ax3 = gca;
VDJidx = find(VDDJtotal>0);
[Rx, Cx] = ind2sub(size(VDDJtotal),VDJidx);
Fx = VDDJtotal(VDJidx);
ScaleFx = MaxDotSize/max(Fx);

scatter(Cx,Rx,Fx*ScaleFx,'fill','MarkerFaceColor','k');
colormap(CustomColor)
set(Ax3,'XTick',[1:length(Xlabels)],'YTick',[1:length(Ylabels)]);
set(Ax3,'XTickLabel',Xlabels,'YTickLabel',Ylabels,'XTickLabelRotation',90);
set(Ax3,'FontName','Arial','FontSize',12)
set(Ax3,'TickLength',[0 0],'box','on')

hold(Ax3,'on')
plot(Ax3,length(Vcount)*[1 1]+1,[0 2*length(Vcount)],'-k','LineWidth',1)
hold(Ax3,'off')
%==========================================================================
%Begin formatting the plots to get the right positions, etc.

%1) Set the x limits
Xlimits = [0 size(VDDJtotal,2)+1];
Ylimits = [0 size(VDDJtotal,1)+1];
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
AllMax = max([VJtotal(:); Dtotal(:)]);
AllMax = AllMax*1.05;
if AllMax <= 0; return; end %****************NOTHING TO PLOT!
set(Ax1,'ylim',[0 AllMax]);
set(Ax2,'xlim',[0 AllMax]);

%Set the X and Y tick labels for frequency, based on maximums
YTick = [0 max(VJtotal(:))];
YLabels = {'0' num2str(max(VJtotal(:)),'%0.0f')};
set(Ax1,'YTick',YTick,'YTickLabel',YLabels);
XTick = [0 max(Dtotal(:))];
XLabels = {'0' num2str(max(Dtotal(:)),'%0.0f')};
set(Ax2,'XTick',XTick,'XTickLabel',XLabels,'XTickLabelRotation',90);

%Adjust the bar colors
Bx1.FaceColor = [0 0 0];
Bx2.FaceColor = [0 0 0];


%resizeSubplots(Gx, 'ScaleVertical', 'y', 'HorzSpacer', 0.01, 'VertSpacer', 0.01, 'FigSpacer', 0.01, ExpPu{:});

FullSaveName = prepSaveTarget(ExpPu{:}, 'SaveSubDir', SaveSubDir, 'SaveExt', '.png', 'DefaultSaveAs', DefaultSaveAs, 'MakeSaveDir', 'y');
savePlot(Gx, ExpPu{:}, 'SaveAs', FullSaveName);

if strcmpi(Visible, 'off')
    close(Gx);
end

varargout{1} = FullSaveName;