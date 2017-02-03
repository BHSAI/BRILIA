
function [H1T, AxT] = plotBarGraph(Edges,Ncount,Normalize,XTickName,LegendName,TitleName)
%Plots the data
FontSize = [10 10 12 10];
if Normalize == 1
    Ncount = Ncount./repmat(sum(Ncount,1),size(Ncount,1),1);
    Yrange = [0 1];
    Yname = 'Norm Freq';
else
    Yrange = [0 max(Ncount)];
    Yname = 'Freq';
end

%Plot the histogram
H1T = bar(Edges(1:end-1),Ncount);
AxT = gca;
set(AxT,'XLim',[Edges(1)-1 Edges(end)]);
set(AxT,'YLim',Yrange);
set(AxT,'Xtick',Edges(1:end-1),'XTickLabelRotation',90); %remmeber to use 1 less tick.
set(AxT,'XTickLabel',XTickName);
modifyPlotLabels(AxT,'',Yname,TitleName,LegendName,FontSize);

%Format the colors
ColorMat = [1 0 0; 0 1 0; 0.2 0.2 1;...
            1 1 0; 0 1 1; 1 0 1];
for j = 1:length(H1T)
    H1T(j).FaceColor = ColorMat(j,:);
end