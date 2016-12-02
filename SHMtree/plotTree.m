[FileName,FilePath] = uigetfile('*.tree','Select the tree file');

%Determine line count first
FID = fopen([FilePath FileName],'r');
NumElem = 0;
while feof(FID) == 0
    TextLine= fgetl(FID);
    if strcmpi(TextLine(1),'>')
        NumElem = NumElem+1;
    end
end

%Create the cell array
frewind(FID)
NewickCell = cell(NumElem,3);
j = 1;
while j <= NumElem
    TextLine= fgetl(FID);
    if strcmpi(TextLine(1),'>')
        NewickCell{j,1} = strrep(TextLine,'>',''); %GroupNum
        NewickCell{j,2} = fgetl(FID); %Newick String
        NewickCell{j,3} = fgetl(FID); %Template count info
        j = j+1;
    end
    
    if feof(FID) == 1
        break
    end
end
fclose(FID);

%Plot the trees
for j = 1:NumElem
%     if isempty(regexp(NewickCell{j,2},'\,','once'))
%         continue
%     end
    TreeObj = phytreeread(NewickCell{j,2}); 
    TreePlot = plot(TreeObj,'Type','square');
    
    %Preformatting the plots
    TreePlot.LeafDots.Marker = 'none';
    TreePlot.BranchDots.Marker = 'none';
    Xlimits = get(TreePlot.axes,'XLim')
    set(TreePlot.axes,'XTick',[-1:1:ceil(Xlimits(2))]);
    AXT = title(NewickCell{j,1});
    set(AXT,'FontSize',14)
        
    %Getting the Leaves, then overlaying scatter plots
    NodeNames = get(TreeObj,'NodeNames');
    NodeCoord = zeros(length(NodeNames),3);
    for k = 1:length(NodeNames)
        NodeCoord(k,1) = TreePlot.BranchLines(k).XData(1);
        NodeCoord(k,2) = TreePlot.BranchLines(k).YData(1);
    end
    
    %Extracting the template count info, and adding to NodeCoord col 3
    TemplateData = regexp(NewickCell{j,3},'\:|\,','split');
    Nrows = length(TemplateData)/2;
    TemplateData = reshape(TemplateData,2,Nrows)'; %{'Names'   'Template Count'}
    for q = 1:size(TemplateData,1)
        TemplateData{q,2} = eval(TemplateData{q,2});
    end
    TemplateLoc = findHeader(NodeNames,TemplateData(:,1));
    DelThis = TemplateLoc == 0;
    TemplateLoc(DelThis) = [];
    TemplateData(DelThis,:) = [];
    NodeCoord(TemplateLoc,3) = cell2mat(TemplateData(:,2));
    KeepThese = NodeCoord(:,3) > 0;
    NodeCoord = NodeCoord(KeepThese,:);
    
    %Overlay the scatterplot
    ScaleFactor = 20;
    hold(TreePlot.axes,'on');
    AX = scatter(TreePlot.axes,NodeCoord(:,1),NodeCoord(:,2),NodeCoord(:,3)*ScaleFactor,[0 0 1],'fill');
    pause
end
% 
% %Extract the key information needed to draw scatter plots
% TreePlot.BranchDots
% TreePlot.branchNodeLabels