%comparing the NT composiiton of the D region.
load('IMGT_Mouse_VDJgene_AllStrain.mat')

NcompF = '';
NcompR = '';
for j = 1:2:size(Dmap,1)
    NcompF = [NcompF Dmap{j,1}];
end
for j = 2:2:size(Dmap,1)
    NcompR = [NcompR Dmap{j,1}];
end

NcompF = '';
for j = 1:1:size(Vmap,1)
    NcompF = [NcompF; Vmap{j,1}(end-2:end)];
    if strcmpi(Vmap{j,1}(end-2:end),'ccc'); pause ; end
end

Vmap(292,:)

N1count = basecount(NcompF);
N2count = basecount(NcompR);
BaseNames = fieldnames(N1count);

N1ct = cell2mat(struct2cell(N1count));
N2ct = cell2mat(struct2cell(N2count));
N1ct = N1ct/sum(N1ct);
N2ct = N2ct/sum(N2ct);


%==========================================================================
Edges = 1:length(BaseNames);
h1 = bar(Edges,[N2ct N1ct],'histc');
%Extract the h1 path width
PatchWidth = diff(unique(h1(1).XData(:,1)));
h1 = bar(Edges-PatchWidth,[N2ct N1ct],'histc'); %Redraw with shift
    h1(1).FaceColor = [1 0 0];
    h1(2).FaceColor = [0 1 0];
    xlim([Edges(1)-PatchWidth Edges(end)+PatchWidth])
    ylim([0 1])
    legend('Forward D','Inverse D');
    ax1 = gca;
    set(ax1,'FontSize',14)
    set(ax1,'TickDir','Out')
    set(ax1,'XTick',Edges)
    set(ax1,'XTickLabel',BaseNames)
    xlabel('# of NTs')
    ylabel('Norm Freq')
    title('N nucleutide insertion composition')
    
set(gcf,'PaperPosition',[0 0 11 8])

OutputFilePre = 'DgeneNTcomp';
saveas(gcf,[OutputFilePre 'CmprDcomp.png']);
save([OutputFilePre 'CmprDcomp.mat'],'N1ct','N2ct');
clear
