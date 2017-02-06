%setPlotTickLabels will extract the plot's current X and Y tick labels, and
%then normalize the decimals in it. Usefule if Matlab yields varying number
%of decimals points.
%
%  setPlotTickLabels(Ax,DecNum) will add DecNum decimal points to current
%  labels.
%  
%  if DecNum = 1
%  if DecNum = [1 2], will add 1 decimal to X, 2decimals to Y 
%  if DecNum = [-1 2], will keep decimal of X, but change only Y. Neg are
%  ignored.

function setPlotTickLabels(Ax,DecNum)
if length(DecNum) == 1
    Xdec = DecNum;
    Ydec = DecNum;
else
    Xdec = DecNum(1);
    Ydec = DecNum(2);
end

if Xdec >= 0
    Xval = get(Ax,'XTick');
    Xstr = cell(1,length(Xval));
    StrForm = ['%0.' num2str(Xdec) 'f'];
    for x = 1:length(Xstr)
        Xstr{x} = sprintf(StrForm,Xval(x));
    end
    set(Ax,'XTickLabel',Xstr);
end

if Ydec >= 0
    Yval = get(Ax,'YTick');
    Ystr = cell(1,length(Yval));
    StrForm = ['%0.' num2str(Ydec) 'f'];
    for y = 1:length(Ystr)
        Ystr{y} = sprintf(StrForm,Yval(y));
    end
    set(Ax,'YTickLabel',Ystr);
end

