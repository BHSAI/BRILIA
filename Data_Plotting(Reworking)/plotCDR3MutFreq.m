function Gx = plotCDR3MutFreq(FreqMat,varargin)
if isempty(varargin)
    AAorder = int2aa(1:20);
else
    AAorder = varargin{1};
    AAmap = aa2int(AAorder);
    
    %Reorder the FreqMat to the New Order
    FreqMat = FreqMat(:,AAmap);
    FreqMat = FreqMat(AAmap,:);
    
end
NTlabel = cell(1,20);
for k = 1:length(NTlabel)
    NTlabel{k} = AAorder(k); 
end
    
FreqMat2 = FreqMat;%./repmat(sum(FreqMat,1),size(FreqMat,1),1);
Crange = [0 max(FreqMat2(:))];

figure
image(FreqMat2,'CDataMapping','scaled')
colormap(jet)
set(gca,'XTick',1:length(AAorder))
set(gca,'XTickLabel',NTlabel)
set(gca,'XAxisLocation','top')
set(gca,'YTick',1:length(AAorder))
set(gca,'YTickLabel',NTlabel)
set(gca,'OuterPosition',[0 0 1 1]);
set(gca,'Clim',Crange);
colorbar

Gx = gcf;