%analyzeN1N2comp will analyze the composition of NTs located in the
%sequences.

%Open file.
[FileName, FilePath] = uigetfile('*.xlsx','Open Excel file');
[~, ~, SampleData] = xlsread([FilePath FileName]);
[SampleData, Header] = filterHeader(SampleData);
NTseq = SampleData(:,3);

%Setup output file names prefix and also decide if to cluster the data
DotLoc = regexp(FileName,'\.');
DotLoc = DotLoc(end);
OutputFilePre = FileName(1:DotLoc(end)-1);

while 1
    ClusterOn = input('Do you want to cluster results? y or n ','s');
    if strcmpi(ClusterOn,'y')
        [GrpNums, GrpNumUnqIdx, ~] = unique(cell2mat(SampleData(:,end)));
        SampleData = SampleData(GrpNumUnqIdx,:);
        break;
    elseif strcmpi(ClusterOn,'n')
        break
    end
end

for k = 1:length(Header)
    if strcmpi(Header{k},'vIndex')
        VidxCol = k;
    elseif strcmpi(Header{k},'n2Index')
        N2idxCol = k;
    elseif strcmpi(Header{k},'dIndex')
        DidxCol = k;
    elseif strcmpi(Header{k},'n1Index')
        N1idxCol = k;
    elseif strcmpi(Header{k},'jIndex')
        JidxCol = k;
    end
end

N1comp = '';
N2comp = '';
for j = 1:size(SampleData,1)
    Seq = SampleData{j,1};
    N2idx = SampleData{j,N2idxCol};
    Didx = SampleData{j,DidxCol};
    N1idx = SampleData{j,N1idxCol};
    Jidx = SampleData{j,JidxCol};
    N1comp = [N1comp Seq(N2idx:Didx-1)];
    N2comp = [N2comp Seq(N1idx:Jidx-1)];
end

N1count = basecount(N1comp);
N2count = basecount(N2comp);
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
    legend('N2','N1');
    ax1 = gca;
    set(ax1,'FontSize',14)
    set(ax1,'TickDir','Out')
    set(ax1,'XTick',Edges)
    set(ax1,'XTickLabel',BaseNames)
    xlabel('# of NTs')
    ylabel('Norm Freq')
    title('N nucleutide insertion composition')
    
set(gcf,'PaperPosition',[0 0 11 8])
saveas(gcf,[OutputFilePre 'CmprN1N2comp.png']);
save([OutputFilePre 'CmprN1N2comp.mat'],'N1ct','N2ct');
clear