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
ClusterInfo1 = cell2mat(VDJdata1(:,[H.SeqNumLoc,H.GrpNumLoc,H.CDR3Loc(2)])); %[SeqNum GrpNum CDR3Len]
ClusterInfo2 = cell2mat(VDJdata2(:,[H.SeqNumLoc,H.GrpNumLoc,H.CDR3Loc(2)])); %[SeqNum GrpNum CDR3Len]

%Calculate the comparison data matrix XYdata
UnqGrpNum = unique(ClusterInfo2(:,2));
XYdata = zeros(length(UnqGrpNum),4);
ClustNumT0 = zeros(size(ClusterInfo1,1),1) > 1;
for y = 1:length(UnqGrpNum)
    
    Grp2Loc = find(UnqGrpNum(y) == ClusterInfo2(:,2));
    Seq2 = ClusterInfo2(Grp2Loc,1);
    [Seq1,Seq1Loc,~] = intersect(ClusterInfo1(:,1),Seq2);
    Grp1Num = ClusterInfo1(Seq1Loc,2);
    UnqGrp1 = unique(Grp1Num);
% 
%     %Refind the full IMGTseqLoc
%     ClustNumT = ClustNumT0;
%     for k = 1:length(UnqGrp1)
%         AddNum = ClusterInfo1(:,2) == UnqGrp1(k);
%         ClustNumT(AddNum) = 1;
%     end
%    Grp1Num = ClusterInfo1(ClustNumT,2);
    
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

%Adding slight perturbations around x and y to better see the dots
Adj = 0.2*(rand([size(XYdata,1) 1])-0.5);
Gx = figure;
Ax = gca;
scatter(Ax,XYdata(:,1) + Adj(:,1),XYdata(:,2)+Adj,30,'ko','fill')
xlim([0 200])
ylim([0 20])
set(get(Ax,'Xlabel'),'string','BRILIA Cluster Size','FontName','Arial','FontSize',18);
set(get(Ax,'Ylabel'),'string','Standard Ovlp Cluster','FontName','Arial','FontSize',18);
set(Ax,'FontName','Arial','FontSize',14);

%Fit a line to the scatter plots
hold(gca,'on');
plot([0 200],[1 1],'k-','LineWidth',1);
hold(gca,'off')

formatPlotSize(Ax,4,4,0);
saveas(gcf,[FilePath, FileNamePre1 '_Tree-ClustDiffIMGT.png']);
close(gcf)

%--------------------------------------------------------------------------
%Plot the correlation in cluster size
Gx = figure;
Ax = gca;
X = XYdata(:,1);
Y = XYdata(:,4)./X;
scatter(X,Y,30,'ko','fill')
xlim([0 200])
ylim([0 1])
set(get(Ax,'Xlabel'),'string','BRILIA Cluster Size','FontName','Arial','FontSize',18);
set(get(Ax,'Ylabel'),'string','Standard Max Ovlp Fr','FontName','Arial','FontSize',18);
set(Ax,'FontName','Arial','FontSize',14);
YtickVal = [0:0.1:1];
YtickLab = cell(length(YtickVal),1);
for k = 1:length(YtickLab)
    YtickLab{k} = num2str(YtickVal(k),'%0.1f');
end
set(Ax,'YTick',YtickVal,'YTickLabel',YtickLab)

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

[Gx,Ax] = formatPlotSize(Ax,4,4,0);
saveas(gca,[FilePath, FileNamePre1 '_Tree-ClustCorrIMGT.png']);

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




%==========================================================================

%==========================================================================

%==========================================================================
%For every unique cluster in IMGT, find how many clusters are in BRILIA
ClusterInfoT = ClusterInfo2;
ClusterInfo2 = ClusterInfo1;
ClusterInfo1 = ClusterInfoT;

%Calculate the comparison data matrix XYdata
UnqGrpNum = unique(ClusterInfo2(:,2));
XYdata = zeros(length(UnqGrpNum),4);
ClustNumT0 = zeros(size(ClusterInfo1,1),1) > 1;
for y = 1:length(UnqGrpNum)
    
    Grp2Loc = find(UnqGrpNum(y) == ClusterInfo2(:,2));
    Seq2 = ClusterInfo2(Grp2Loc,1);
    [Seq1,Seq1Loc,~] = intersect(ClusterInfo1(:,1),Seq2);
    Grp1Num = ClusterInfo1(Seq1Loc,2);
    UnqGrp1 = unique(Grp1Num);
% 
%     %Refind the full IMGTseqLoc
%     ClustNumT = ClustNumT0;
%     for k = 1:length(UnqGrp1)
%         AddNum = ClusterInfo1(:,2) == UnqGrp1(k);
%         ClustNumT(AddNum) = 1;
%     end
%     Grp1Num = ClusterInfo1(ClustNumT,2);
%     
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

%Adding slight perturbations around x and y to better see the dots
Adj = 0.2*(rand([size(XYdata,1) 1])-0.5);
Gx = figure;
Ax = gca;
scatter(Ax,XYdata(:,1) + Adj(:,1),XYdata(:,2)+Adj,30,'ko','fill')
xlim([0 200])
ylim([0 20])
set(get(Ax,'Xlabel'),'string','Standard Cluster Size','FontName','Arial','FontSize',18);
set(get(Ax,'Ylabel'),'string','BRILIA Ovlp Cluster','FontName','Arial','FontSize',18);
set(Ax,'FontName','Arial','FontSize',14);

%Fit a line to the scatter plots
hold(gca,'on');
plot([0 200],[1 1],'k-','LineWidth',1);
hold(gca,'off')

formatPlotSize(Ax,4,4,0);
saveas(gcf,[FilePath, FileNamePre1 '_Tree-ClustDiffBRILIA.png']);
close(gcf)

%--------------------------------------------------------------------------
%Plot the correlation in cluster size
Gx = figure;
Ax = gca;
X = XYdata(:,1);
Y = XYdata(:,4)./X;
scatter(X,Y,30,'ko','fill')
xlim([0 200])
ylim([0 1])
set(get(Ax,'Xlabel'),'string','Standard Cluster Size','FontName','Arial','FontSize',18);
set(get(Ax,'Ylabel'),'string','BRILIA Max Ovlp Fr','FontName','Arial','FontSize',18);
set(Ax,'FontName','Arial','FontSize',14);
YtickVal = [0:0.1:1];
YtickLab = cell(length(YtickVal),1);
for k = 1:length(YtickLab)
    YtickLab{k} = num2str(YtickVal(k),'%0.1f');
end
set(Ax,'YTick',YtickVal,'YTickLabel',YtickLab)

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

[Gx,Ax] = formatPlotSize(Ax,4,4,0);
saveas(gca,[FilePath, FileNamePre1 '_Tree-ClustCorrBRILIA.png']);














% 
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
%Comparing Vmut Freq
H.VmutLoc = findHeader(VDJheader,'vMutCt_Germline');
CDR3length = VDJdata1{1,H.CDR3Loc(2)};
VframeLength = length(VDJdata1{1,H.SeqLoc}) - CDR3length;

%Get the V frame and CDR3 mut counts, with repsect to germline
for w = 1:2
    if w == 1
        VDJdataT = VDJdata1; 
    else
        VDJdataT = VDJdata2;
    end
    MutData = zeros(size(VDJdata1,1),2); %[Vmut and CDR3mut]
    GrpNum = cell2mat(VDJdataT(:,H.GrpNumLoc));
    UnqGrp = unique(GrpNum);
    for y = 1:length(UnqGrp)
        IdxLoc = find(GrpNum == UnqGrp(y));
        RefSeq = VDJdataT{IdxLoc(1),H.RefSeqLoc};
        RefVframe = RefSeq(1:VframeLength);
        RefCDR3 = RefSeq(VframeLength+1:VframeLength+CDR3length);

        for k = 1:length(IdxLoc)
            CurSeq = VDJdataT{IdxLoc(k),H.SeqLoc};
            CurVframe = CurSeq(1:VframeLength);
            CurCDR3 = CurSeq(VframeLength+1:VframeLength+CDR3length);

            MutData(IdxLoc(k),1) = sum(RefVframe ~= CurVframe);
            MutData(IdxLoc(k),2) = sum(RefCDR3 ~= CurCDR3);
        end
    end
    
    if w == 1
        MutData1 = MutData;
    else
        MutData2 = MutData;
    end
end

EdgesVmut = [0:0.03:0.21]; %[1:10]/VframeLength; %Determine bin sizes based on length of 104C framework / 10, or CDR3 nt length /10
Vmut1 = MutData1(:,1)/VframeLength; %cell2mat(VDJdata1(:,H.VmutLoc));%./cell2mat(VDJdata1(:,H.LengthLoc(1)));
Vmut2 = MutData2(:,1)/VframeLength; %cell2mat(VDJdata2(:,H.VmutLoc));%./cell2mat(VDJdata2(:,H.LengthLoc(1)));
BinsVmut1 = histc(Vmut1,EdgesVmut);
BinsVmut2 = histc(Vmut2,EdgesVmut);
b6 = bar(EdgesVmut,[BinsVmut1,BinsVmut2],'histc');
b6(1).FaceColor = [0 0 1];
b6(2).FaceColor = [1 0 0];
Xticks = EdgesVmut; %[-0.04:0.04:EdgesVmut(end)+0.04];
%xlim([Xticks(1) Xticks(end)])
xlim([0 0.21])
ylim([0 80])
set(gca,'XTick',Xticks)
set(gca,'YTick',[0:20:80])
set(get(gca,'xlabel'),'string','Vframe Mut Fr','FontSize',14,'FontName','Arial')
set(get(gca,'ylabel'),'string','Frequency','FontSize',14,'FontName','Arial')
set(gca,'FontName','Arial','FontSize',16);
formatPlotSize(gca,5,2); 
saveas(gcf,[FilePath, FileNamePre1 '_Vframemut.png']);

EdgesCDR3mut = [0:0.03:0.21]; %[1:10]/CDR3length; %Determine bin sizes based on length of 104C framework / 10, or CDR3 nt length /10
CDR3mut1 = MutData1(:,2)/CDR3length; %cell2mat(VDJdata1(:,H.VmutLoc));%./cell2mat(VDJdata1(:,H.LengthLoc(1)));
CDR3mut2 = MutData2(:,2)/CDR3length; %cell2mat(VDJdata2(:,H.VmutLoc));%./cell2mat(VDJdata2(:,H.LengthLoc(1)));
BinsCDR3mut1 = histc(CDR3mut1,EdgesCDR3mut);
BinsCDR3mut2 = histc(CDR3mut2,EdgesCDR3mut);
b6 = bar(EdgesCDR3mut,[BinsCDR3mut1,BinsCDR3mut2],'histc');
b6(1).FaceColor = [0 0 1];
b6(2).FaceColor = [1 0 0];
Xticks = EdgesCDR3mut; %[-0.04:0.04:EdgesVmut(end)+0.04];
%xlim([Xticks(1) Xticks(end)])
xlim([0 0.21])
ylim([0 80])
set(gca,'XTick',Xticks)
set(gca,'YTick',[0:20:80])
set(get(gca,'xlabel'),'string','CDR3 Mut Fr','FontSize',14,'FontName','Arial')
set(get(gca,'ylabel'),'string','Frequency','FontSize',14,'FontName','Arial')
set(gca,'FontName','Arial','FontSize',16);
formatPlotSize(gca,5,2); 
saveas(gcf,[FilePath, FileNamePre1 '_CDR3mut.png']);



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

