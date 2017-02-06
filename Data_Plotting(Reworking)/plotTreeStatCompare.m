%This will open IMGT and BRILIA annotation, and then compare the extra
%subcluster vs BRILIA cluster size.

[VDJdata1,VDJheader,FileName1,FilePath] = openSeqData; %IMGT
[VDJdata2,VDJheader,FileName2,FilePath] = openSeqData; %BRILIA
DotLoc = find(FileName1 == '.');
FileNamePre1 = FileName1(1:DotLoc(end)-1);
DotLoc = find(FileName2 == '.');
FileNamePre2 = FileName2(1:DotLoc(end)-1);

H = getHeaderVar(VDJheader);

%Get tree stats
ClusterInfoA = cell2mat(VDJdata1(:,[H.SeqNumLoc,H.GrpNumLoc,H.CDR3Loc(2)])); %[SeqNum GrpNum CDR3Len]
ClusterInfoB = cell2mat(VDJdata2(:,[H.SeqNumLoc,H.GrpNumLoc,H.CDR3Loc(2)])); %[SeqNum GrpNum CDR3Len]

for w = 1:2
    if w == 1
        ClusterInfo1 = ClusterInfoA;
        ClusterInfo2 = ClusterInfoB;
        XaxisLabA = 'Cluster Size';
        YaxisLabA = 'Standard Ovlp Cluster';
        TitleLabA = 'BRILIA Clusters';
        XaxisLabB = 'Cluster Size';
        YaxisLabB = 'Standard Corresponding Cluster Size';
        TitleLabB = 'BRILIA Clusters';
    
    else
        ClusterInfo1 = ClusterInfoB;
        ClusterInfo2 = ClusterInfoA;
        XaxisLabA = 'Cluster Size';
        YaxisLabA = 'BRILIA Ovlp Cluster';
        TitleLabA = 'Standard Clusters';
        XaxisLabB = 'Cluster Size';
        YaxisLabB = 'BRILIA Corresponding Cluster Size';
        TitleLabB = 'Standard Clusters';
    end

    %Calculate the comparison data matrix XYdata
    UnqGrpNum = unique(ClusterInfo2(:,2));
    XYdata = zeros(length(UnqGrpNum),4);
    ClustNumT0 = zeros(size(ClusterInfo1,1),1) > 1;
    for y = 1:length(UnqGrpNum)
        %Find the Grp1Clusters that overlap with Grp2Clusters
        Grp2Loc = find(UnqGrpNum(y) == ClusterInfo2(:,2));
        Seq2 = ClusterInfo2(Grp2Loc,1);
        [Seq1,Seq1Loc,~] = intersect(ClusterInfo1(:,1),Seq2);
        Grp1Num = ClusterInfo1(Seq1Loc,2);
        UnqGrp1 = unique(Grp1Num);

        %Refind the full IMGTseqLoc
        ClustNumT = ClustNumT0;
        for k = 1:length(UnqGrp1)
            AddNum = ClusterInfo1(:,2) == UnqGrp1(k);
            ClustNumT(AddNum) = 1;
        end
        Grp1Num = ClusterInfo1(ClustNumT,2);

        MaxCt = 0;
        for k = 1:length(UnqGrp1)
            CurCt = sum(UnqGrp1(k) == Grp1Num);
            if CurCt > MaxCt
                MaxCt = CurCt;
            end
        end

        %Find the maximum size overlap
        XYdata(y,1) = length(Grp2Loc); %unq seq2 per grp2
        XYdata(y,2) = length(UnqGrp1); %unq cluster ct per grp2
        XYdata(y,3) = ClusterInfo2(Grp2Loc(1),3); %CDR3 length
        XYdata(y,4) = MaxCt; %unq1 seq in largest subset grp
    end

    if w == 1
        %PULL OUT the very divergent ones, with multiple sub clusters
        KeepThese = find(XYdata(:,2) == max(XYdata(:,2)));
        XYdataT = XYdata(KeepThese,:);
        XYdiff = XYdataT(:,1) - XYdataT(:,4);
        MaxDiffLoc = KeepThese(find(XYdiff == max(XYdiff)));
        MaxDiffLoc = MaxDiffLoc(1); %GrpNumber of interest

        GrpLoc = find(MaxDiffLoc == ClusterInfo2(:,2));
        Seq2 = ClusterInfo2(GrpLoc,1);
        [Seq1,Seq1Loc,~] = intersect(ClusterInfo1(:,1),Seq2);
        ClustNum = ClusterInfo1(Seq1Loc,2);
        UnqGrp1 = unique(ClustNum); %IMGT's grp loc

        ExtractIMGTLoc = [];
        for k = 1:length(UnqGrp1)
            AddLoc = find(UnqGrp1(k) == ClusterInfo1(:,2));
            ExtractIMGTLoc = [ExtractIMGTLoc; AddLoc];
        end
        ExtractMavricLoc = GrpLoc;

        SaveData1 = [VDJheader; VDJdata1(ExtractIMGTLoc,:)];
        SaveData2 = [VDJheader; VDJdata2(ExtractMavricLoc,:)];
        xlswrite('WorstDiff_IMGT.xlsx',SaveData1);
        xlswrite('WorstDiff_Mavric.xlsx',SaveData2);
        
        %------------------------------------------------------------------
        %Save dot locations for this worst case example onto the two othe plots
        GrpSave{2} = unique(cell2mat(VDJdata1(ExtractIMGTLoc,H.GrpNumLoc))); %The numberfor GrpSave is flipped, due to the way define X and Y
        GrpSave{1} = unique(cell2mat(VDJdata2(ExtractMavricLoc,H.GrpNumLoc)));
    end
    
    %Plot the example dot on top of the plots
    SelGrp = GrpSave{w};
    SelLoc = zeros(size(XYdata,1),1)>1;
    for z = 1:length(SelGrp)
        SelLoc = (UnqGrpNum == SelGrp(z)) | SelLoc;
    end
    
    %Adding slight perturbations around x and y to better see the dots
    Adj = 0;%0.2*(rand([size(XYdata,1) 1])-0.5);
    
    %----------------------------------------------------------------------
    %Plotting Cluster2 Size on X, Cluster1 Ovlp Ct on Y
    
    %Setting the variables and labels
    X = XYdata(:,1) + Adj;
    Y = XYdata(:,2) + Adj;
    Xticks = [0:40:200];
    Yticks = [0:2:20];
    DotClr = zeros(size(X,1),3);
    DotClr(SelLoc,1) = 1;
    
    %Reoder the dots so that the red ones are always on top
    [~,SortIdx] = sort(SelLoc);
    X = X(SortIdx);
    Y = Y(SortIdx);
    DotClr = DotClr(SortIdx,:);
    
        %Set the xy ticks, limits, labels
        Gx1 = figure;
        Ax1 = gca;
        scatter(Ax1,X,Y,40,DotClr,'o','fill')

        Xlims = [Xticks(1) Xticks(end)];
        Ylims = [Yticks(1) Yticks(end)];
        Xlabs = cell(1,length(Xticks));
        Ylabs = cell(1,length(Yticks));
        for p = 1:length(Xticks)
            Xlabs{p} = num2str(Xticks(p),'%0.0f');
        end
        for p = 1:length(Yticks)
            Ylabs{p} = num2str(Yticks(p),'%0.0f');
        end
        set(Ax1,'Xlim',Xlims,'Ylim',Ylims,'XTick',Xticks,'YTick',Yticks,'XTickLabel',Xlabs,'YTickLabel',Ylabs);        

        %Set the X, Y, Title texts
        set(get(Ax1,'Xlabel'),'string',XaxisLabA,'FontName','Arial','FontSize',14);
        set(get(Ax1,'Ylabel'),'string',YaxisLabA,'FontName','Arial','FontSize',14);
        set(get(Ax1,'Title'),'string',TitleLabA,'FontName','Arial','FontSize',14);
        set(Ax1,'FontName','Arial','FontSize',14,'Box','on');

    %Fit a line to the scatter plots
    hold(Ax1,'on');
    plot([0 200],[1 1],'k-','LineWidth',1);
    hold(Ax1,'off')
    formatPlotSize(Ax1,4,4,0);

    %----------------------------------------------------------------------
    %Plotting Cluster2 Size on X, Cluster1 Corresponding Cluster Y

    %Setting the variables and labels
    X = XYdata(:,1) + Adj;
    Y = XYdata(:,4) + Adj;
    Xticks = [0:40:200];
    Yticks = [0:40:200];
    DotClr = zeros(size(X,1),3);
    DotClr(SelLoc,1) = 1;
    
    %Reoder the dots so that the red ones are always on top
    [~,SortIdx] = sort(SelLoc);
    X = X(SortIdx);
    Y = Y(SortIdx);
    DotClr = DotClr(SortIdx,:);

        %Set the xy ticks, limits, labels
        Gx2 = figure;
        Ax2 = gca;
        scatter(Ax2,X,Y,40,DotClr,'o','fill')

        Xlims = [Xticks(1) Xticks(end)];
        Ylims = [Yticks(1) Yticks(end)];
        Xlabs = cell(1,length(Xticks));
        Ylabs = cell(1,length(Yticks));
        for p = 1:length(Xticks)
            Xlabs{p} = num2str(Xticks(p),'%0.0f');
        end
        for p = 1:length(Yticks)
            Ylabs{p} = num2str(Yticks(p),'%0.0f');
        end
        set(Ax2,'Xlim',Xlims,'Ylim',Ylims,'XTick',Xticks,'YTick',Yticks,'XTickLabel',Xlabs,'YTickLabel',Ylabs);        

        %Set the X, Y, Title texts
        set(get(Ax2,'Xlabel'),'string',XaxisLabB,'FontName','Arial','FontSize',14);
        set(get(Ax2,'Ylabel'),'string',YaxisLabB,'FontName','Arial','FontSize',14);
        set(get(Ax2,'Title'),'string',TitleLabB,'FontName','Arial','FontSize',14);
        set(Ax2,'FontName','Arial','FontSize',14,'Box','on');
        
    %Fit a line to the scatter plots
    hold(gca,'on');
    plot([0 200],[0 200],'k--','LineWidth',1);
    hold(gca,'off')
    formatPlotSize(Ax2,4,4,0);
    
    %Save the plots
    if w == 1
        saveas(Gx1,[FilePath, FileNamePre1 '_ClustDiff.png']);
        saveas(Gx2,[FilePath, FileNamePre1 '_ClustCorr.png']);
    else
        saveas(Gx1,[FilePath, FileNamePre2 '_ClustDiff.png']);
        saveas(Gx2,[FilePath, FileNamePre2 '_ClustCorr.png']);
    end
    close(Gx1)
    close(Gx2)
end

% %Fit a line to the scatter plots
% Pfits = polyfit(X,Y,1);
% Yfits = polyval(Pfits,X);
% Yresid = Y - Yfits;
% SSresid = sum(Yresid.^2);
% SStotal = (length(Y)-1) * var(Y);
% Rsq = 1 - SSresid/SStotal;
% hold(gca,'on');
% plot(X,Yfits,'k','LineWidth',1,'Color',[1 0 0 ]);
% plot([0 200],[0 200],'-','LineWidth',1,'Color',[0.5 0.5 0.5]);
% hold(gca,'off')
% text(100,180,sprintf('y = %0.2fx+%0.2f (R^2 = %0.2f)',Pfits(1),Pfits(2),Rsq),'HorizontalAlignment','center','FontName','Arial','FontSize',14,'Color',[1 0 0])

% % 
% % 
% % %==========================================================================
% % [FileName1,FilePath1] = uigetfile('*.mat','Select IMGT mat for TreeStats')
% % [FileName2,FilePath2] = uigetfile('*.mat','Select BRILIA mat for TreeStats')
% % load([FilePath1,FileName1]);
% % Tree1 = TreeStats;
% % load([FilePath2,FileName2]);
% % Tree2 = TreeStats;
% 

% Tree1 = findTreeStats(VDJdata1,VDJheader);
% Tree2 = findTreeStats(VDJdata2,VDJheader);

% %Comparing MaxNodeLevel
% EdgesNode = [0:1:15];
% BinsNode1 = histc(Tree1(:,2),EdgesNode);
% BinsNode2 = histc(Tree2(:,2),EdgesNode);
% b1 = bar(EdgesNode,[BinsNode1,BinsNode2]);
% b1(1).FaceColor = [0 0 1];
% b1(2).FaceColor = [1 0 0];
% xlim([0 EdgesNode(end)])
% set(gca,'XTick',EdgesNode)
% set(get(gca,'xlabel'),'string','MaxNodeLevel','FontSize',14,'FontName','Arial')
% set(get(gca,'ylabel'),'string','Frequency','FontSize',14,'FontName','Arial')
% formatPlotSize(gca,5,5); 
% saveas(gcf,[FilePath, FileNamePre1 '_Tree-MaxNode.png']);
% 
% %Comparing MaxChildPerNode
% EdgesChild = [0:1:max([Tree1(:,3) ; Tree2(:,3)])];
% BinsChild1 = histc(Tree1(:,3),EdgesChild);
% BinsChild2 = histc(Tree2(:,3),EdgesChild);
% b2 = bar(EdgesChild,[BinsChild1,BinsChild2]);
% b2(1).FaceColor = [0 0 1];
% b2(2).FaceColor = [1 0 0];
% xlim([0 30])
% set(gca,'XTick',EdgesChild)
% set(get(gca,'xlabel'),'string','MaxChilderPerNode','FontSize',14,'FontName','Arial')
% set(get(gca,'ylabel'),'string','Frequency','FontSize',14,'FontName','Arial')
% formatPlotSize(gca,5,5); 
% saveas(gcf,[FilePath, FileNamePre1 '_Tree-MaxChildPNode.png']);
% 
% %Comparing MaxExtentOfMutaiton
% EdgesMut = [0:0.02:max([Tree1(:,4) ; Tree2(:,4)])];
% BinsMut1 = histc(Tree1(:,4),EdgesMut);
% BinsMut2 = histc(Tree2(:,4),EdgesMut);
% b3 = bar(EdgesMut,[BinsMut1,BinsMut2]);
% b3(1).FaceColor = [0 0 1];
% b3(2).FaceColor = [1 0 0];
% xlim([-0.02 EdgesMut(end)])
% set(gca,'XTick',EdgesMut)
% set(get(gca,'xlabel'),'string','MaxMutationExtent','FontSize',14,'FontName','Arial')
% set(get(gca,'ylabel'),'string','Frequency','FontSize',14,'FontName','Arial')
% formatPlotSize(gca,5,5); 
% saveas(gcf,[FilePath, FileNamePre1 '_Tree-MaxMutExtent.png']);
% 
% %Comparing NumberUnqSeq
% EdgesSeq = [0:1:max([Tree1(:,5) ; Tree2(:,5)])];
% BinsSeq1 = histc(Tree1(:,5),EdgesSeq);
% BinsSeq2 = histc(Tree2(:,5),EdgesSeq);
% b4 = bar(EdgesSeq,[BinsSeq1,BinsSeq2]);
% b4(1).FaceColor = [0 0 1];
% b4(2).FaceColor = [1 0 0];
% %xlim([0 EdgesSeq(end)])
% xlim([0 30])
% set(gca,'XTick',EdgesSeq)
% set(get(gca,'xlabel'),'string','UnqSeq','FontSize',14,'FontName','Arial')
% set(get(gca,'ylabel'),'string','Frequency','FontSize',14,'FontName','Arial')
% formatPlotSize(gca,5,5); 
% saveas(gcf,[FilePath, FileNamePre1 '_Tree-UnqSeq.png']);
% 
% %Comparing MaxHamDistExplored
% EdgesDist = [0:1:max([Tree1(:,6) ; Tree2(:,6)])];
% BinsDist1 = histc(Tree1(:,6),EdgesDist);
% BinsDist2 = histc(Tree2(:,6),EdgesDist);
% b5 = bar(EdgesDist,[BinsDist1,BinsDist2]);
% b5(1).FaceColor = [0 0 1];
% b5(2).FaceColor = [1 0 0];
% xlim([0 EdgesDist(end)])
% set(gca,'XTick',EdgesDist)
% set(get(gca,'xlabel'),'string','MaxHamDist','FontSize',14,'FontName','Arial')
% set(get(gca,'ylabel'),'string','Frequency','FontSize',14,'FontName','Arial')
% formatPlotSize(gca,5,5); 
% saveas(gcf,[FilePath, FileNamePre1 '_Tree-MaxHamDist.png']);
% 
% %Comparing NumberNodeperGroup
% EdgesNodeCt = [0:1:max([Tree1(:,7) ; Tree2(:,7)])];
% BinsNodeCt1 = histc(Tree1(:,7),EdgesNodeCt);
% BinsNodeCt2 = histc(Tree2(:,7),EdgesNodeCt);
% b0 = bar(EdgesNodeCt,[BinsNodeCt1,BinsNodeCt2]);
% b0(1).FaceColor = [0 0 1];
% b0(2).FaceColor = [1 0 0];
% xlim([0 15])
% set(gca,'XTick',EdgesNodeCt)
% set(get(gca,'xlabel'),'string','NodeCt','FontSize',14,'FontName','Arial')
% set(get(gca,'ylabel'),'string','Frequency','FontSize',14,'FontName','Arial')
% formatPlotSize(gca,5,5); 
% saveas(gcf,[FilePath, FileNamePre1 '_Tree-NodeCt.png']);
% 
% 
% %==========================================================================
% %Comparing annotation level differences in mutations
% 

EdgesVmut = [0:20]; %Determine bin sizes based on length of 104C framework / 10, or CDR3 nt length /10
VDJmut1 = sum(cell2mat(VDJdata1(:,H.VmutLoc:H.JmutLoc)),2);%./cell2mat(VDJdata1(:,H.LengthLoc(1)));
VDJmut2 = sum(cell2mat(VDJdata2(:,H.VmutLoc:H.JmutLoc)),2);%./cell2mat(VDJdata2(:,H.LengthLoc(1)));
BinsVmut1 = histc(VDJmut1,EdgesVmut);
BinsVmut2 = histc(VDJmut2,EdgesVmut);
b6 = bar(EdgesVmut,[BinsVmut1,BinsVmut2],'histc');
b6(1).FaceColor = [0 0 1];
b6(2).FaceColor = [1 0 0];

Xlabs = cell(size(EdgesVmut));
for j = 1:length(EdgesVmut)
    Xlabs{j} = num2str(EdgesVmut(j),'%0.0f');
end
Xticks = EdgesVmut+b6(1).XData(3,1); %You add 0.5 to center the labels
xlim([EdgesVmut(1) 15])
%ylim([0 120])
set(gca,'XTick',Xticks,'XTickLabel',Xlabs)
set(get(gca,'xlabel'),'string','Germline-Child Seq Mutation','FontSize',14,'FontName','Arial')
set(get(gca,'ylabel'),'string','Frequency','FontSize',14,'FontName','Arial')
formatPlotSize(gca,4,4); 
saveas(gcf,[FilePath, FileNamePre1 '_VDJmut.png']);



% 
% EdgesVmut = [0:20]; %Determine bin sizes based on length of 104C framework / 10, or CDR3 nt length /10
% Vmut1 = cell2mat(VDJdata1(:,H.VmutLoc));%./cell2mat(VDJdata1(:,H.LengthLoc(1)));
% Vmut2 = cell2mat(VDJdata2(:,H.VmutLoc));%./cell2mat(VDJdata2(:,H.LengthLoc(1)));
% BinsVmut1 = histc(Vmut1,EdgesVmut);
% BinsVmut2 = histc(Vmut2,EdgesVmut);
% b6 = bar(EdgesVmut,[BinsVmut1,BinsVmut2],'histc');
% b6(1).FaceColor = [0 0 1];
% b6(2).FaceColor = [1 0 0];
% Xticks = EdgesVmut; %[-0.04:0.04:EdgesVmut(end)+0.04];
% xlim([Xticks(1) Xticks(end)])
% ylim([0 120])
% set(gca,'XTick',Xticks)
% set(get(gca,'xlabel'),'string','Vmut','FontSize',14,'FontName','Arial')
% set(get(gca,'ylabel'),'string','Frequency','FontSize',14,'FontName','Arial')
% formatPlotSize(gca,5,3); 
% saveas(gcf,[FilePath, FileNamePre1 '_Vmut.png']);
% 
% %Comparing Dmut Freq
% H.DmutLoc = findHeader(VDJheader,'dMutCt_Germline');
% EdgesDmut = [0:20];%[0:1/8:0.3];
% Dmut1 = cell2mat(VDJdata1(:,H.DmutLoc));%./cell2mat(VDJdata1(:,H.LengthLoc(3)));
% Dmut2 = cell2mat(VDJdata2(:,H.DmutLoc));%./cell2mat(VDJdata2(:,H.LengthLoc(3)));
% BinsDmut1 = histc(Dmut1,EdgesDmut);
% BinsDmut2 = histc(Dmut2,EdgesDmut);
% b7 = bar(EdgesDmut,[BinsDmut1,BinsDmut2],'histc');
% b7(1).FaceColor = [0 0 1];
% b7(2).FaceColor = [1 0 0];
% Xticks = EdgesDmut; %[-0.04:0.04:EdgesDmut(end)+0.04];
% xlim([Xticks(1) Xticks(end)])
% ylim([0 120])
% 
% set(gca,'XTick',EdgesDmut)
% set(get(gca,'xlabel'),'string','Dmut','FontSize',14,'FontName','Arial')
% set(get(gca,'ylabel'),'string','Frequency','FontSize',14,'FontName','Arial')
% formatPlotSize(gca,5,3); 
% saveas(gcf,[FilePath, FileNamePre1 '_Dmut.png']);
% 
% 
% %Comparing Jmut Freq
% H.JmutLoc = findHeader(VDJheader,'jMutCt_Germline');
% EdgesJmut = [0:20]; %;[0:1/16:0.3];
% Jmut1 = cell2mat(VDJdata1(:,H.JmutLoc));%./cell2mat(VDJdata1(:,H.LengthLoc(5)));
% Jmut2 = cell2mat(VDJdata2(:,H.JmutLoc));%./cell2mat(VDJdata2(:,H.LengthLoc(5)));
% BinsJmut1 = histc(Jmut1,EdgesJmut);
% BinsJmut2 = histc(Jmut2,EdgesJmut);
% b8 = bar(EdgesJmut,[BinsJmut1,BinsJmut2],'histc');
% b8(1).FaceColor = [0 0 1];
% b8(2).FaceColor = [1 0 0];
% Xticks = EdgesJmut; %[-0.04:0.04:EdgesJmut(end)+0.04];
% xlim([Xticks(1) Xticks(end)])
% ylim([0 120])
% 
% set(gca,'XTick',EdgesJmut)
% set(get(gca,'xlabel'),'string','Jmut','FontSize',14,'FontName','Arial')
% set(get(gca,'ylabel'),'string','Frequency','FontSize',14,'FontName','Arial')
% formatPlotSize(gca,5,3); 
% saveas(gcf,[FilePath, FileNamePre1 '_Jmut.png']);

% %==========================================================================
% %Compare the CDR3 mutation frequencies now, IMGT vs BRILIA, using maximum
% %deviation.
% 
% Output1 = findVvsCDR3mut(VDJdata1,VDJheader);
% [Gx3,Ax3,PDFmat3] = plotVvsCDR3mut(Output1);
% 
% Output2 = findVvsCDR3mut(VDJdata2,VDJheader);
% [Gx4,Ax4,PDFmat4] = plotVvsCDR3mut(Output2);
% 
% %Set the scale the same, by setting PDFmat(1,1) to be same.
% MaxVal = max([PDFmat3(2:end) PDFmat4(2:end)])+2;
% PDFmat3(1,1) = MaxVal;
% PDFmat4(1,1) = MaxVal;
% Ix3 = get(Ax3,'children');
% Ix4 = get(Ax4,'children');
% set(Ix3,'CData',PDFmat3);
% set(Ix4,'CData',PDFmat4);
% 
% %Set the colorbar info
% c3 = colorbar('peer',Ax3,'east','FontSize',12);
% ColorBarLim = get(c3,'Limits');
% ColorBarTicks = [0:ColorBarLim(2)/10:ColorBarLim(2)];
% ColorBarSpacing = ColorBarTicks/ColorBarLim(2)*max(PDFmat3(:));
% ColorBarTickLabel = cell(size(ColorBarSpacing));
% for j = 1:length(ColorBarSpacing)
%     ColorBarTickLabel{j} = sprintf('%0.0f',ColorBarSpacing(j));
% end
% set(c3,'Ticks',ColorBarTicks,'TickLabels',ColorBarTickLabel);
% 
% %Set the colorbar info
% c4 = colorbar('peer',Ax4,'east','FontSize',12);
% ColorBarLim = get(c4,'Limits');
% ColorBarTicks = [0:ColorBarLim(2)/10:ColorBarLim(2)];
% ColorBarSpacing = ColorBarTicks/ColorBarLim(2)*max(PDFmat4(:));
% ColorBarTickLabel = cell(size(ColorBarSpacing));
% for j = 1:length(ColorBarSpacing)
%     ColorBarTickLabel{j} = sprintf('%0.0f',ColorBarSpacing(j));
% end
% set(c4,'Ticks',ColorBarTicks,'TickLabels',ColorBarTickLabel);
% 
% saveas(Gx3,[FilePath, FileNamePre1 '_VvsCDR3mut.png']);
% saveas(Gx4,[FilePath, FileNamePre2 '_VvsCDR3mut.png']);

