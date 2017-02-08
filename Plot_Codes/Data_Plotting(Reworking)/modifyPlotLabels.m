%modifyPlotLabels will modify a series of X or Y axis labels, Title,
%Legends. Note that this will only overwrite. Empty values are skipped.
%This could also be used to modify fontsize only, by setting all other
%parameters fixed.

%FontSize = [10 10 14 10] will apply font size to X,Y,Title,Legend,
%respectively.
function modifyPlotLabels(Ax,Xname,Yname,TitleName,LegendName,FontSize)
if isempty(FontSize)
    FontSize = [10 10 10 10];
end

if ~isempty(Xname)
    xlabel(Ax,Xname);
end

if ~isempty(Yname)
    ylabel(Ax,Yname);
end

if ~isempty(TitleName)
    title(Ax,TitleName);
end

if ~isempty(LegendName)
    H4 = legend(Ax,LegendName);
    set(H4,'FontSize',FontSize(4)); %Since legend creates new obj, this is here so you adjust font.
else
    H4 = findobj(gcf,'tag','legend');
    set(H4,'FontSize',FontSize(4));
end

set(get(Ax,'Xlabel'),'FontSize',FontSize(1));
set(get(Ax,'Ylabel'),'FontSize',FontSize(2));
set(get(Ax,'Title'),'FontSize',FontSize(3));

