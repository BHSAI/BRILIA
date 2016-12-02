%masterPlot is a script organizing all the data anlysis to be done.
%
% All scripts are arranged to do the following
% 1) Extract the data from the xlsx files using the "find____" functions
% 2) Plot the extracted data using the "plot____" functions
% 3) Format the plots for publishable quality "format____" functions

NumFiles = input('How many files are you comparing? ');
FileName = cell(NumFiles,1);
FilePath = cell(NumFiles,1);
FileData = cell(NumFiles,2);
for j = 1:NumFiles
    [FileName{j}, FilePath{j}] = uigetfile('*.xlsx',sprintf('Open File #%d',j));
    [~, ~, SampleData] = xlsread([FilePath{j}, FileName{j}]);
    [FileData{j,1}, FileData{j,2}] = filterHeader(SampleData);
end

[OutputFile, OutputPath] = uiputfile('*.*','Where do you want to save output files?');
ClusterThis = 0;
while 1
    ClusterOn = input('Do you want to cluster results? y or n ','s');
    if strcmpi(ClusterOn,'y')
        OutputFile = [OutputFile '.Clustered'];
        ClusterThis = 1;
        break;
    elseif strcmpi(ClusterOn,'n')
        break
    end
end

%Create the File Header Legend
FileLegend = cell(1,size(FileName,1));
for j = 1:length(FileName);
    FileLegend{j} = sprintf('File%02d',j);
end

%==========================================================================
%Comparing deletions and insertions

AllHeaderEval = {'vDeletion', 'n2AlignLength', 'd5Deletion', 'd3Deletion', 'n1AlignLength', 'jDeletion'};%, 'dMatchCt'};

for k = 1:length(AllHeaderEval)
    HeaderEval = AllHeaderEval(k);
    ColData = cell(1,length(FileName));
    for j = 1:length(FileName)
        ColData(j) = findGeneMods(FileData{j,1},FileData{j,2},HeaderEval,ClusterThis);
    end
    [~, Ax1] = plotGeneMods(ColData,[0:25],'','norm');
    if length(FileName) == 1
        TitleName = sprintf('%s, Navg = %0.2f, Nstd = %0.2f',HeaderEval{1},mean(ColData{1}),std(ColData{1}));
        modifyPlotLabels(Ax1,'','',TitleName,FileLegend,[10 10 14 10]); %Use HeaderEval as Title
    else
        modifyPlotLabels(Ax1,'','',HeaderEval,FileLegend,[10 10 14 10]); %Use HeaderEval as Title
    end
    
    %Print figure as png
    set(gcf,'PaperPosition',[0 0 5 5]);
    HeaderNospace = HeaderEval{1};
    HeaderNospace(HeaderNospace == ' ') = [];
    saveas(gcf,[OutputPath OutputFile '_' HeaderNospace '.png']);
end
% 
% %==========================================================================
% %Comparing VDJ family usage, inidividually
% 
% FamData = cell(1,length(FileName));
% for j = 1:length(FileName)
%     FamData(:,j) = findGeneFamilyUsage(FileData{j,1},FileData{j,2},ClusterThis);
% end
% [~, Ax2] = plotGeneFamilyUsage(FamData,FileLegend,'norm'); %This spits out a full 1x3 plot
% 
% %Print figure as png
% Pre = 'VDJ';
% for j = 1:length(Ax2)
%     FigNum = get(Ax2{j},'Parent');
%     set(FigNum,'PaperPosition',[0 0 5 5]);
%     Post = sprintf('_FamilyUsage%s.png',Pre(j));
%     saveas(FigNum,[OutputPath OutputFile Post]);
% end
% %--------------------------------------------------------------------------
% %Comparing VDJ family usage, as V-D-J combo scatter plots
% 
% Ax3 = plotGeneFamilyCombo(FamData,FileLegend,3000);
% 
% %Print figure as png
% for j = 1:length(FileName)
%     FigNum = get(Ax3{j,1},'Parent');
%     FileLabel = sprintf('File%0d',j);
%     set(FigNum,'PaperPosition',[0 0 10 10]);
%     saveas(FigNum,[OutputPath OutputFile '_FamilyCombo' FileLabel '.png']);
% end
% 
% %==========================================================================
% %Comparing CDR3 lengths
% 
% HeaderEval = {'cdr3Length'};
% ColData = cell(1,length(FileName));
% for j = 1:length(FileName)
%     ColData(j) = findGeneMods(FileData{j,1},FileData{j,2},HeaderEval,ClusterThis);
%     ColData{j} = ColData{j}/3;
% end
% [~, Ax1] = plotGeneMods(ColData,[0:25],'');
% modifyPlotLabels(Ax1,'','',HeaderEval,FileLegend,[10 10 14 10]);  
% 
% %Print figure as png
% set(gcf,'PaperPosition',[0 0 5 5]);
% HeaderNospace = HeaderEval{1};
% HeaderNospace(HeaderNospace == ' ') = [];
% saveas(gcf,[OutputPath OutputFile '_' HeaderNospace '.png']);