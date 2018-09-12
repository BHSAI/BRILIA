%getTrunkLenDist will get the trunk length distributions. 

function getTreeLenDist()
FileList = openFileDialog('.csv', 'Open csv files', 'multiselect', 'on');
for f = 1:length(FileList)
	getTreeLenDist_Local(FileList{f});
end

function getTreeLenDist_Local(InFileName)

[VDJdata, VDJheader, FileName, FilePath] = openSeqData(InFileName);
[H, ~, Chain] = getAllHeaderVar(VDJheader);

if strcmpi(Chain, 'none')
    for j = [3:10 12:length(VDJheader)]
        VDJheader{j} = ['H-' VDJheader{j}];
    end
    H = getAllHeaderVar(VDJheader);
end

GrpNum = cell2mat(VDJdata(:, H.GrpNumLoc));
UnqGrpNum = unique(GrpNum);

LenStats = zeros(length(UnqGrpNum), 7); %GrpCt, CDR3len, SeqLen, MinMiss, MinTrunkLen%, MaxMiss, MaxTrunkLen%
for j = 1:length(UnqGrpNum)
    Idx = find(UnqGrpNum(j) == GrpNum);
    RefSeq = VDJdata{Idx(1), H.RefSeqLoc};
    CDR3Len = VDJdata{Idx(1), H.CDR3Loc(2)};
    if ischar(CDR3Len)
        CDR3Len = str2double(CDR3Len);
    end
    MaxMissCt = -1;
    for k = 1:length(Idx)
        Seq = VDJdata{Idx(k), H.SeqLoc};
        MissCt = sum(~((Seq == RefSeq) | Seq == 'N' | RefSeq == 'N'));
        if k == 1
            MinMissCt = MissCt;
        end
        if MissCt > MaxMissCt
            MaxMissCt = MissCt;
        end
    end
    LenStats(j, :) = [length(Idx) CDR3Len  length(RefSeq)  MinMissCt  MinMissCt/length(RefSeq) MaxMissCt MaxMissCt/length(RefSeq)];
end

%For CDF

%Remove singletones
LenStats(LenStats(:, 1) == 1, :) = [];

%Lowerbound cdf
MinSHM = sort(LenStats(:,5));
MinY = [1:length(MinSHM)]/length(MinSHM);

%Upperbound cdf
MaxSHM = sort(LenStats(:,7));
MaxY = [1:length(MaxSHM)]/length(MaxSHM);

plot(MinSHM, MinY, 'b', MaxSHM, MaxY, 'r', 'LineWidth', 1);

DotLoc = find(FileName == '.');
title(strrep(FileName(1:DotLoc(1)-1), '_', '\_'));

xlabel('% mutation');
ylabel('clonotype frequency');
resizeSubplots(gca);
% 
% OutFileName = fullfile(FilePath, [FileName(1:DotLoc(end)-1), 'TreeLenDist.png']);
% savePlot(gcf, 'SaveAs', OutFileName, 'Format', 'png');
% % 
% 
% 
% %For histograms
% Bins = 0.00:1/125:0.15; 
% 
% Nmin = histc(LenStats(:, end-2), Bins);
% Nmax = histc(LenStats(:, end), Bins);
% X = repelem(Bins, [1 2*ones(1, numel(Bins)-2) 3]);
% Ymin = repelem(Nmin, 2*ones(1, numel(Bins)));
% Ymax = repelem(Nmax, 2*ones(1, numel(Bins)));
% Lx = plot(X, Ymin, 'b', X, Ymax, 'r');
% Lx(1).LineWidth = 2;
% Lx(2).LineWidth = 2;
% 
% DotLoc = find(FileName == '.');
% 
% 
% title(strrep(FileName(1:DotLoc(1)-1), '_', '\_'));
% 
% BinTxt = cellfun(@(x)sprintf('%0.1f', x*100), num2cell(Bins), 'UniformOutput', false);
% 
% set(gca, 'XTickLabel', BinTxt(1:3:end), 'XTick', Bins(1:3:end), 'fontsize', 16, 'linewidth', 1);
% xlabel('% mutation');
% ylabel('clonotype frequency');
% resizeSubplots(gca);
xlim([0 0.3]);
OutFileName = fullfile(FilePath, [FileName(1:DotLoc(end)-1), 'TreeLenDist.png']);
savePlot(gcf, 'SaveAs', OutFileName, 'Format', 'png');

