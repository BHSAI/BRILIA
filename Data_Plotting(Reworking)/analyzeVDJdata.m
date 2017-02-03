%analyzeVDJresults will take an output EXCEL file and then plot data
[FileName, FilePath] = uigetfile('*.xlsx;*.mat','Select VDJ results file');
[~, ~, FileExt] = fileparts([FilePath FileName]);
if FileExt == '.xlsx'
    [~, ~, TableData] = xlsread([FilePath FileName]);
    HeaderT = TableData(1,:);
    TableData(1,:) = [];
elseif FileExt == '.mat'
    load([FilePath FileName],'Results');
    Results = buildVDJclassifier(Results);
    Results = cleanVDJresults(Results);
    Results = buildVDJclassifier(Results);
    TableData = convertResults2Table(Results);
end

%Find the filename pre
DotLoc = find(FileName == '.');
DotLoc = DotLoc(end);
FileNamePre = FileName(1:DotLoc-1);

%Find CDR3 AA and NT lengths. Append and save.
if  size(TableData,2) == 26
    CDR3data = cell(size(TableData,1),2);
    for k = 1:size(TableData,1)
        CDR3aa = findCDR3(TableData{k,2},TableData{k,3});
        CDR3data{k,1} = CDR3aa;
        CDR3data{k,2} = length(CDR3aa)*3;
    end
    TableData = cat(2,TableData,CDR3data);
    HeaderT = cat(2,HeaderT,{'CDR3' 'CDR3 len'});
    xlswrite([FilePath FileName],cat(1,HeaderT,TableData));
end

%Open the VDJ Gene database in .mat file, if it exist. 
DefaultFile = 'IMGT_Mouse_VDJgene_AllStrain.mat';
if exist(DefaultFile,'file') > 0
    load(DefaultFile);
else
    error('Need "IMGT_Mouse_VDJgene_C57BL6.mat". Look at generateIMGTmat.m code.');
end

%Convert all family name string to matrix
FamCols = [12 17 22];
for j = 1:size(TableData,1)
    for k = 1:3
        if ischar(TableData{j,FamCols(k)});
            TableData{j,FamCols(k)} = eval(TableData{j,FamCols(k)});
        end
    end
end

%Filter out non-funcitonal ones
BadOnes = (cell2mat(TableData(:,end)) == 0);
TableData(BadOnes,:) = [];

%Make Histogram of individual family usage
GrpNums = cell2mat(TableData(:,6));
[~,GrpUnqLoc,~] = unique(GrpNums);
Vbin = zeros(size(Vmap,1),1);
Dbin = zeros(size(Dmap,1),1);
Jbin = zeros(size(Jmap,1),1);
for j = 1:length(GrpUnqLoc)
    VbinLoc = TableData{GrpUnqLoc(j),FamCols(1)};
    DbinLoc = TableData{GrpUnqLoc(j),FamCols(2)};
    JbinLoc = TableData{GrpUnqLoc(j),FamCols(3)};
    Vbin(VbinLoc) = Vbin(VbinLoc) + 1/length(VbinLoc);
    Dbin(DbinLoc) = Dbin(DbinLoc) + 1/length(DbinLoc);
    Jbin(JbinLoc) = Jbin(JbinLoc) + 1/length(JbinLoc);
end

%FIGURE 1
figure(1);

%Group by gene family name (not gene level)
GeneFamNum = Vmap(:,5);
[VmapUnq,VmapUnqIdx,UnqLoc] = unique(GeneFamNum);
Vbin2 = zeros(size(VmapUnq));
for k = 1:max(UnqLoc)
    Vloc = find(UnqLoc == k);
    Vbin2(k) = sum(Vbin(Vloc));    
end
subplot(2,2,1)
bar(Vbin2);
Gx1 = gca;
set(Gx1,'XTickLabel',VmapUnq);
set(Gx1,'XTickLabelRotation',90);
set(Gx1,'XTick',[1:length(VmapUnq)]);

GeneFamNum = Jmap(:,5);
[JmapUnq,JmapUnqIdx,UnqLoc] = unique(GeneFamNum);
Jbin2 = zeros(size(JmapUnq));
for k = 1:max(UnqLoc)
    Vloc = find(UnqLoc == k);
    Jbin2(k) = sum(Jbin(Vloc));    
end
subplot(2,2,3)
bar(Jbin2);
Gx3 = gca;
set(Gx3,'XTickLabel',JmapUnq);
set(Gx3,'XTickLabelRotation',90);
set(Gx3,'XTick',[1:length(JmapUnq)]);

GeneFamNum = Dmap(:,5);
[DmapUnq,DmapUnqIdx,UnqLoc] = unique(GeneFamNum);
Dbin2 = zeros(size(DmapUnq));
for k = 1:max(UnqLoc)
    Vloc = find(UnqLoc == k);
    Dbin2(k) = sum(Dbin(Vloc));    
end
DbinFwd = Dbin2(1:2:end);
DbinRev = Dbin2(2:2:end);
DmapUnqFwd = DmapUnq(1:2:end);
DmapUnqRev = DmapUnq(2:2:end);

subplot(2,2,2)
bar(DbinFwd);
Gx2 = gca;
set(Gx2,'XTickLabel',DmapUnqFwd);
set(Gx2,'XTickLabelRotation',90);
set(Gx2,'XTick',[1:length(DmapUnqFwd)]);

subplot(2,2,4)
bar(DbinRev);
Gx4 = gca;
set(Gx4,'XTickLabel',DmapUnqRev);
set(Gx4,'XTickLabelRotation',90);
set(Gx4,'XTick',[1:length(DmapUnqRev)]);
set(Gx4,'YLim',get(Gx2,'YLim'))

saveas(gcf,[FileNamePre '.GeneFam.png']);


%FIGURE 2

%Filter to get only unique groups
TableData = TableData(GrpUnqLoc,:);

figure(2)
Edges = [0:15];
%Vdeletion
ColNum = 16;
N = histc(cell2mat(TableData(:,ColNum)),Edges);
subplot(2,2,1);
bar(Edges,N,'histc');
title('Vgene 3'' Deletion')

%Jdeletion
ColNum = 24;
N = histc(cell2mat(TableData(:,ColNum)),Edges);
subplot(2,2,3);
bar(Edges,N,'histc');
title('Jgene 5'' Deletion')

%5Ddeletion
ColNum = 19;
N = histc(cell2mat(TableData(:,ColNum)),Edges);
subplot(2,2,2);
bar(Edges,N,'histc');
title('Dgene 5'' Deletion')

%3Ddeletion
ColNum = 21;
N = histc(cell2mat(TableData(:,ColNum)),Edges);
subplot(2,2,4);
bar(Edges,N,'histc');
title('Dgene 3'' Deletion')

%FIGURE 3
figure(3)
XlocV = [1:size(Vbin2)];
YlocJ = [1:size(Jbin2)];
ZlocD = [1:size(Dbin2)];

saveas(gcf,[FileNamePre '.GeneDel.png']);

close all

%Make 3D scatter plot of V D J genes
GrpNums = cell2mat(TableData(:,6));
[~,GrpUnqLoc,~] = unique(GrpNums);
[VmapUnq,~,VmapIdx] = unique(Vmap(:,5));

Dbin = zeros(size(Dmap,1),1);
[DmapUnq,~,DmapIdx] = unique(Dmap(:,5));

Jbin = zeros(size(Jmap,1),1);
[JmapUnq,~,JmapIdx] = unique(Jmap(:,5));

VDJbin = zeros(length(VmapUnq),length(DmapUnq),length(JmapUnq));
for j = 1:length(GrpUnqLoc)
    VbinLoc = VmapIdx(TableData{GrpUnqLoc(j),FamCols(1)});
    DbinLoc = DmapIdx(TableData{GrpUnqLoc(j),FamCols(2)});
    JbinLoc = JmapIdx(TableData{GrpUnqLoc(j),FamCols(3)});
    VDJbin(VbinLoc,DbinLoc,JbinLoc) = VDJbin(VbinLoc,DbinLoc,JbinLoc) + 1/(length(VbinLoc)*length(DbinLoc)*length(JbinLoc));
end
BinIdx = find(VDJbin > 0);
[X,Y,Z] = ind2sub(size(VDJbin),BinIdx);

ScaleFactor = 100;
BinFreq = VDJbin(BinIdx)*ScaleFactor;

ColorMat = rand(size(BinFreq,1),3);
scatter3(X,Y,Z,BinFreq,ColorMat,'filled')

for Jnum = 1:4
    figure(Jnum)
    scatter(X(Z == Jnum),Y(Z == Jnum),BinFreq(Z == Jnum),'filled')
    xlim([0 16])
    ylim([0 12])
    AX1 = gca;
    set(AX1,'XTickLabel',VmapUnq);
    set(AX1,'XTick',[1:length(VmapUnq)]);
    set(AX1,'XTickLabelRotation',90);
    set(AX1,'YTickLabel',DmapUnq);
    set(AX1,'YTick',[1:length(DmapUnq)]);
    title(sprintf('J%01d Family',Jnum));
    saveas(gcf,sprintf('%s.Scatter VDJ_J%1d.png',FileNamePre,Jnum));
end


%Planetary Plot! Shows VDJ family combo + CDR3 avg length with STD
%variation as a circular ring. 


