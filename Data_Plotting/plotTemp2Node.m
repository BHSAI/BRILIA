%plotTemp2Node plots the group template count on y axis, and the number of
%unique sequences on x axis. the linearity shows something interesting...

[FileNames,FilePath] = uigetfile('*.xlsx','Open VDJdata files','multiselect','on');
if ischar(FileNames)
    FileNames = {FileNames};
end
FileNames = sort(FileNames);

%For multiple bar graph comparison
Nedges = [0:1:30];
Bins = zeros(length(Nedges),length(FileNames));

for f = 1:length(FileNames)
    FileName = FileNames{f};
    DotLoc = find(FileName == '.');
    FileNamePre = FileName(1:DotLoc(end)-1);

    [VDJdata,NewHeader] = openSeqData([FilePath FileNames{f}]);
    

    getHeaderVar
    
    %Renumber group number, assuming groups are already grouped together in
    %the file.
    GrpN = 1;
    SetNum = VDJdata{1,GrpNumLoc};
    for j = 1:size(VDJdata,1)
       CurNum = VDJdata{j,GrpNumLoc};
       if CurNum ~= SetNum
           GrpN = GrpN + 1;
           SetNum = CurNum;
       end
       VDJdata{j,GrpNumLoc} = GrpN;
    end
    
    %Extract the group data = [GrpNum GrpSize SumTemp MaxTemp Vnum Dnum Jnum]
    GrpNum = cell2mat(VDJdata(:,GrpNumLoc));
    GrpNumUnq = unique(GrpNum);
    GrpTempMat = zeros(length(GrpNumUnq),4);
    
    for y = 1:length(GrpNumUnq)
        GrpIdx = find(GrpNum == GrpNumUnq(y));
        TempCt = cell2mat(VDJdata(GrpIdx,TemplateLoc));
        GrpTempMat(y,1) = GrpNumUnq(y);
        GrpTempMat(y,2) = length(GrpIdx);
        GrpTempMat(y,3) = sum(TempCt);
        GrpTempMat(y,4) = max(TempCt);
        
        for k = 1:3
            Xname = VDJdata{GrpIdx,FamLoc(k)};
            ColonLoc = find(Xname == ':');
            if ~isempty(ColonLoc)
                Xname = Xname(ColonLoc+1:end);
            end
            Xname = strrep(Xname,' ','');
            BarLoc = find(Xname =='|');
            if ~isempty(BarLoc)
                Xname = Xname(1:BarLoc-1);
            end
            if strcmpi(Xname(1),'r')
                Xnum = Xname(6:7);
                Dir = -1;
            else
                Xnum = Xname(5:6);
                Dir = 1;
            end
            DigLoc = regexp(Xnum,'\d');
            Xnum = eval(Xnum(DigLoc))*Dir;
            
            GrpTempMat(y,4+k) = Xnum;
        end
    end
    
    %Getting the count data
    Bins(:,f) = histc(GrpTempMat(:,2),Nedges);
    
    
    %Create the color matrix fo the V,D,J
    
    Gx = figure();    
    for d = 1:max(GrpTempMat(:,5))
        subplot(4,4,d)
        Ax = gca;
        hold(Ax,'on')
        PlotFormat = {'x','o','s'};
        PlotIdx = GrpTempMat(:,5) == d;
        scatter(Ax,GrpTempMat(PlotIdx,2),GrpTempMat(PlotIdx,4),10,'fill');%,GrpTempMat(PlotIdx,5)
        colormap('jet')
        
        TitleName = strrep(FileNames{f},'_','\_');
        title(TitleName)
        xlabel('Num Unq Seq in Grp')
        ylabel('Max Templat Ct in Grp')
        xlim([0 400])
        ylim([0 400])     
        hold(Ax,'off')
        title(sprintf('Vfamily %2d',d))
    end
    set(Gx,'PaperPosition',[0 0 12 12])
    saveas(Gx,[FilePath FileNamePre 'V.png']);

    Gx = figure();
    b = 1;
    for d = min(GrpTempMat(:,6)):max(GrpTempMat(:,6))
        if d == 0; continue; end
        subplot(4,3,b)
        Ax = gca;
        hold(Ax,'on')
        PlotFormat = {'x','o','s'};
        PlotIdx = GrpTempMat(:,6) == d;
        scatter(Ax,GrpTempMat(PlotIdx,2),GrpTempMat(PlotIdx,4),10,'fill');%,GrpTempMat(PlotIdx,5)
        colormap('jet')

        TitleName = strrep(FileNames{f},'_','\_');
        title(TitleName)
        xlabel('Num Unq Seq in Grp')
        ylabel('Max Templat Ct in Grp')
        xlim([0 400])
        ylim([0 400])     
        hold(Ax,'off')
        title(sprintf('dfamily %2d',d))
        b = b+1;
    end
    set(Gx,'PaperPosition',[0 0 12 9])
    saveas(Gx,[FilePath FileNamePre 'D.png']);
    
    Gx = figure();
    for d = 1:max(GrpTempMat(:,7))
        subplot(1,4,d)
        Ax = gca;
        hold(Ax,'on')
        PlotFormat = {'x','o','s'};
        PlotIdx = GrpTempMat(:,7) == d;
        scatter(Ax,GrpTempMat(PlotIdx,2),GrpTempMat(PlotIdx,4),10,'fill');%,GrpTempMat(PlotIdx,5)
        colormap('jet')

        TitleName = strrep(FileNames{f},'_','\_');
        title(TitleName)
        xlabel('Num Unq Seq in Grp')
        ylabel('Max Templat Ct in Grp')
        xlim([0 400])
        ylim([0 400])     
        hold(Ax,'off')
        title(sprintf('jfamily %2d',d))
    end
    set(Gx,'PaperPosition',[0 0 12 3])
    saveas(Gx,[FilePath FileNamePre 'J.png']);

% 
%     %Identify the top picks in either MaxCount or UnqSeqCount
%     TopPicks = 4;
%     if size(GrpTempMat,1) < TopPicks
%         TopPicks = size(GrpTempMat,1);
%     end
%     TempCt = cell2mat(VDJdata(:,TemplateLoc));
%     
%     %Find the top unique sequence counts
%     GrpTempMat = sortrows(GrpTempMat,-2); %UnqSeqCount
%     Top2UnqSeq = GrpTempMat(1:TopPicks,:);
%     Top2UnqSeqCDR3 = cell(TopPicks,1);
%     for q = 1:TopPicks
%         IdxLoc1 = Top2UnqSeq(q,1) == GrpNum;       
%         IdxLoc2 = TempCt == max(TempCt(IdxLoc1));
%         IdxLoc = find(IdxLoc1 & IdxLoc2);
%         IdxLoc = IdxLoc(1);
%         CDR3seq = VDJdata{IdxLoc,CDR3Loc(1)};
%         if ~ischar(CDR3seq)
%             CDR3seq = 'noCDR3';
%         end
%         Top2UnqSeqCDR3{q} = CDR3seq;
%     end
%     
%     %Find the top max template counts
%     GrpTempMat = sortrows(GrpTempMat,-4);
%     Top2TempCt = GrpTempMat(1:TopPicks,:);
%     Top2TempCtCDR3 = cell(TopPicks,1);
%     for q = 1:TopPicks
%         IdxLoc1 = Top2TempCt(q,1) == GrpNum;       
%         IdxLoc2 = TempCt == max(TempCt(IdxLoc1));
%         IdxLoc = find(IdxLoc1 & IdxLoc2);
%         IdxLoc = IdxLoc(1);
%         CDR3seq = VDJdata{IdxLoc,CDR3Loc(1)};
%         if ~ischar(CDR3seq)
%             CDR3seq = 'noCDR3';
%         end
%         Top2TempCtCDR3{q} = CDR3seq;
%     end
%     
%     for k = 1:TopPicks
%         if f == 1
%             RotAngle1 = 90;
%             RotAngle2 = 90;
%         else
%             RotAngle1 = atand(Top2UnqSeq(k,4)/Top2UnqSeq(k,2));
%             RotAngle2 = atand(Top2TempCt(k,4)/Top2TempCt(k,2));
%         end
%         text(Top2UnqSeq(k,2),Top2UnqSeq(k,4),Top2UnqSeqCDR3{k},'Rotation',RotAngle1);
%         text(Top2TempCt(k,2),Top2TempCt(k,4),Top2TempCtCDR3{k},'Rotation',RotAngle2);        
%     end
%     
%     FileName = FileNames{f};
%     DotLoc = find(FileName,'.');
%     formatPlotSize(gca,5,5);
%     SaveName = [FileName(1:DotLoc(end)-1) '.png'];
%     saveas(gcf,[FilePath SaveName])    
end


