%script for performing analysis on ebVLP data donMultVarTest and
%donDiffTest for differences
%
%  For cmprRep, having file names with a group identifying will group-wise comparison of repertoires.
%  Ex: specifying "G1 Grp1 G2 Grp2" will mean group 1 will have the name "Grp1" in the file name.
%                  G1                           G2                          G3      
%                  Grp1_SeqFile1.BRILIAvN.csv   Grp2_SeqFile1.BRILIAvN.csv  Grp3_SeqFile1.BRILIAvN.csv
%                  Grp1_SeqFile2.BRILIAvN.csv   Grp2_SeqFile2.BRILIAvN.csv  Grp3_SeqFile2.BRILIAvN.csv

function runGroupAnalysis(varargin)
% P = inputParser;
% addOptional(P, 'FileNames', '', @(x) isempty(x) || iscell(x));
% addParameter(P, 'Alpha', 0.05, @(x) x >= 0 && x < 1);
% addParameter(P, 'Tail', 'B', @(x) startsWith(upper(x), {'L', 'R', 'B'}));
% addParameter(P, 'TestName', 'Groupwise', @ischar);
% parse(P, varargin{:});
% Alpha = P.Results.Alpha;
% Tail = P.Results.Tail;
% TestName = P.Results.TestName;
% 
% VaccineDir = findVaccineDir;
% EbVLP = getBriliaFiles(fullfile(VaccineDir, 'EbolaVLP*Spleen/*'), 'multiselect');
% VEEV = getBriliaFiles(fullfile(VaccineDir, 'VEEV*Spleen/*'), 'multiselect');
% 
% O = Datastore(VEEV, 'GroupBy', 'dir');

[~, O] = getVaccineFiles(Option);

O.getDiversity;
O.saveData('Diversity')

% O.getClonotypeData;
% O.plotFreqData;
% O.savePlot;

% 
% 
% 
% %O.getConvergence
% 
% % O = Datastore(VEEV, 'GroupBy', 'file', 'GroupFormat', 'Grp#');
% 
% 
% Group = O.Group;
% FileNames = O.Files;
% 
% FileNameOnly = cell(size(FileNames));
% for f = 1:length(FileNames)
%     [~, FileNameOnly{f}] = parseFileName(FileNames{f});
% end
% FileDir = fileparts(fileparts(FileNames{1}));
% 
% S = openMultSeqData(FileNames);
% S = renameSameGene(S, 'mouse'); %rename degenerate gene names
% 
% % printSummary
% 
% %Count the frequencies of clonotypes for various properties
% RepFilt = {'AC', 'SC', 'BC'};
% PropFilt = {'H-CDR3_Length', 'H-V_GeneName', 'H-J_GeneName'};
% for s = 1:length(RepFilt)
%     S = countDRF(S, RepFilt{s}, 0);
%     S = countSHMfreq(S, RepFilt{s});
%     S = countRepProp(S, RepFilt{s}, PropFilt{:});
% end
% % S = calcRepDiversity(S);
% % 
% % %Plot all the diveristy stuff
% % SaveDivDir = fullfile(FileDir, 'RepDiversity', filesep);
% % if ~isdir(SaveDivDir)
% %     mkdir(SaveDivDir);
% % end
% % for s = 1:length(S)
% %     plotConfetti(S(s).Diversity.Confetti);
% %     DivData = [fieldnames(S(s).Diversity) struct2cell(S(s).Diversity)];
% %     TxtCell = cell(1, size(DivData, 1)-1);
% %     for a = 1:length(TxtCell)
% %         TxtCell{a} = sprintf('%s = %0.3f ', DivData{a, :});
% %     end
% %     Ax = get(gca);
% %     Tx = text(Ax.XLim(2), Ax.YLim(2), TxtCell, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontName', 'Arial', 'FontSize', 8);
% %     Tx.BackgroundColor = [0 0 0];
% %     Tx.Color = [1 1 1];
% %     SaveName = strrep(FileNameOnly{s}, '.BRILIAv3.csv', '.Confetti.png');
% %     savePlot(gcf, 'SaveAs', fullfile(SaveDivDir, SaveName))
% % end
% 
% %Get each field names storing frequency data
% Fields = fieldnames(S);
% UndLoc = cellfun(@(x) contains(x, '_'), Fields);
% Fields = unique(cellfun(@(x) x(find(x == '_', 1, 'first')+1:end), Fields(UndLoc), 'unif', false));
% 
% %INTER-repertoire comparison between groups and conditions
% FilesToCombine = cell(1, 100);
% ff = 1;
% 
% for r = 1:length(RepFilt)
%     for f = 1:length(Fields)
%         FieldName = [RepFilt{r} '_' Fields{f}];
%         if contains(FieldName, 'SC_SHM'); continue; end %useless data
% 
%         %Conform the dist label and freq w/ respect to Naive all clonotype AC
%         C = {S.(FieldName)};
%         if contains(Fields{f}, 'SHM')
%             C = conformDist(C{:});
%             C = rebinData(C, [0 0.0001 0.0001:0.01:0.1]);
%             C(2, :) = [];
%             C{1, 1} = 0;
%             XDEC = 3;
%             YDEC = 2;
%         elseif contains(Fields{f}, 'CDR3_Length')
%             continue
%             C(CtrlLoc) = {S(CtrlLoc).(['AC_' Fields{f}])};
%             C = conformDist(C{:});
%             C = rebinData(C, 4:21);
%             C(:, 1) = num2cell(4:20);
%             XDEC = 0;
%             YDEC = 0;
%         else
%             C(CtrlLoc) = {S(CtrlLoc).(['AC_' Fields{f}])};
%             C = conformDist(C{:});
%             XDEC = -1;
%             YDEC =  0;
%         end
% 
%         BootSize = sum(cell2mat(C(:, 2:end)), 1);
% 
%         if ~contains(Fields{f}, 'SHM')
%             [C, ~, Scalor] = scaleMultDist(C, Group, CtrlGrp);
%         else %Normalized to total, since you want fractional values %of of population data!!!!
%             Scalor = ones(1, size(C, 2) - 1);
%             for k = 2:size(C, 2)
%                 C(:, k) = num2cell(cell2mat(C(:, k)) / sum(cell2mat(C(:, k))));
%             end
%         end
% 
%         %Prepare data for nonparametric groupwise test
%         VarName = C(:, 1);
%         Data = cell2mat(C(:, 2:end));
%         [Pall, Gall, MED, CI] = groupwiseTest(Data, Group, 'Alpha', Alpha, 'Tail', Tail, 'Iter', 200, 'Size', BootSize);
%         PairH = Gall(:, 3:end)' <= Alpha;
% 
%         %Show top 50 
%         if size(Data, 1) > 50
%             DataSum = max(Data, [], 2);
%             [~, Idx] = sort(DataSum, 'descend');
%             ShowIdx = sort(Idx(1:50));
%         else
%             ShowIdx = [];
%         end
% 
%         MaxData = applyall('max', Data, Group); 
%         MinData = applyall('min', Data, Group);
%         MidData = (MaxData + MinData) / 2;
%         DistData = MaxData - MidData;
%         MidData(isnan(MidData)) = 0;
%         DistData(isnan(DistData)) = 0;
% 
%         %Make the plot
%         plotErrorBox(MidData, 'STD', DistData, 'PairH', PairH, 'XTickLabel', VarName, 'IndivData', Data, 'Group', Group, 'ShowIdx', ShowIdx);
%         TitleName = sprintf('%ss_%s', RepFilt{r}, strrep(Fields{f}, '_', ' '));
%         title(strrep(TitleName, '_', '\_'));
%         FigH = 6;
%         FigW = 0.25*size(MidData, 1);
%         if FigW <  4; FigW =  4; end
%         if FigW > 15; FigW = 15; end
%         resizeFigure(gcf, FigW, FigH);
%         setPlotTickDecimal(gca, XDEC, YDEC);
%         if contains(Fields{f}, 'SHM')
%             XTickLabel = get(gca, 'XTickLabel');
%             XTickLabel{1} = 'No SHM';
%             set(gca, 'XTickLabel', XTickLabel);
%         end
%         drawnow;
%         resizeSubplots;
% 
%         %Save the image and the data
%         SaveNamePre = strrep(TitleName, ' ', '_');
%         SaveDir = fullfile(FileDir, TestName);
%         if ~exist(SaveDir, 'dir')
%             mkdir(SaveDir);
%         end
%         savePlot(gcf, 'SaveAs', fullfile(SaveDir, [SaveNamePre '.png']), 'dpi', 300);
% 
% %             UG = num2cell(unique(Grouping));
% %             MidHeader = cellfun(@(x) sprintf('G%d_Middle', x), UG, 'unif', false);
% %             StdHeader = cellfun(@(x) sprintf('G%d_Range', x), UG, 'unif', false);
% %             
% %             PairPHeader = cell(1, size(Pall, 1));
% %             for k = 1:length(PairPHeader)
% %                 PairPHeader{k} = sprintf('Subj Diff p-val: %d v %d', Pall(k, 1), Pall(k, 2)); 
% %             end
% %             PairP = Pall(:, 3:end)';
% % 
%         GroupHeader = cell(1, length(Group));
%         for k = 1:length(GroupHeader)
%             GroupHeader{k} = sprintf('%d : %s', Group(k), FileNameOnly{k});
%         end
% 
%         GroupPHeader = cell(1, size(Gall, 1));
%         for k = 1:length(GroupPHeader)
%             GroupPHeader{k} = sprintf('%d v %d p-val ', Gall(k, 1), Gall(k, 2)); 
%         end
%         GroupP = Gall(:, 3:end)';
% 
%         TableHeader = horzcat(Fields{f}, GroupHeader, GroupPHeader);
%         TableScalor = ['Scalor' num2cell(Scalor) repmat({' '}, 1, length(TableHeader) - length(Scalor) - 1)];
%         TableData = [VarName num2cell([Data GroupP])];
%         writeDlmFile([TableHeader; TableScalor; TableData], fullfile(SaveDir, [SaveNamePre '.csv']), ',');
% 
%         FilesToCombine{ff} = fullfile(SaveDir, [SaveNamePre '.csv']);
%         ff = incr(ff);
%     end
%     close all
% end
% 
% combineDlmFile(FilesToCombine(1:ff-1), fullfile(SaveDir, 'AllFreq.csv'));
% delete(FilesToCombine{1:ff-1})
