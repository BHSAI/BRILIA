%extractGroupNumbers will look at the excel files and return the number of
%groups found per file. Good for comparing the hamming distance vs group
%number data.


FileList = dir('*.xlsx*');
GroupCount = cell(length(FileList),2);
for j = 1:size(FileList)
    FileName = FileList(j).name;
    [~,~,SampleData] = xlsread(FileName);
    [SampleData, Header] = filterHeader(SampleData);
    GroupLoc = findHeader(Header,'GroupNum');
    GroupMat = cell2mat(SampleData(:,GroupLoc));
    
    GroupCount{j,1} = FileName;
    if isempty(GroupMat(1)); continue; end
    GroupCount{j,2} = max(GroupMat);
end


GroupCount = sortrows(GroupCount,-2);
GroupNum = cell2mat(GroupCount(:,2));
plot([1:length(GroupNum)],GroupNum,'ok');
xlabel('HammingDistance');
ylabel('Number of Clusters');
set(gca,'FontSize',12,'XTick',[1:length(GroupNum)])

set(gcf,'PaperPosition',[0 0 5 5]);

[OutputFile, OutputPath] = uiputfile('*.png','Save Hamming Distance plot as?');
saveas(gcf,[OutputPath OutputFile '_HammingDis.png']);
