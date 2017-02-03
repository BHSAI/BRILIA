function Gx = plotHotSpots2(Amat,Cmat,Gmat,Tmat)

%To simplify the plotting, rearrange matrices
HotMat = {Amat{1} Cmat{1} Gmat{1} Tmat{1}};
AllMat = {Amat{2} Cmat{2} Gmat{2} Tmat{2}};

%Plotting all hot spots
for q = 1:2
    Gx{q} = figure;
    for k = 1:4
        subplot(1,4,k)
        if q == 1
            CurMat = HotMat{k};
        else
            CurMat = AllMat{k};
        end
        MidLoc = ceil(size(CurMat,2)/2);
        CurMat = CurMat / max(CurMat(:,MidLoc));
            

        %Calculate the position for the texts
        TextXpos = [1:size(CurMat,2)];
        TextYpos = cumsum(CurMat,1) - CurMat/2;

        H1 = bar(CurMat','stack');
        H1(1).FaceColor = [1.0 0.2 0.2];
        H1(2).FaceColor = [0.2 1.0 0.2];
        H1(3).FaceColor = [0.4 0.4 1.0];
        H1(4).FaceColor = [0.7 0.7 0.7];
% 
%         ylim([0,1.1])
        
        %Making the labels
        Labels = cell(1,size(CurMat,2));
        for y = 1:length(Labels)
            Labels{y} = num2str(y-MidLoc);
        end
        
        set(gca,'TickLength',[0 0])
        set(gca,'XTickLabel',Labels)
        set(gca,'OuterPosition',[0.25*(k-1) 0 0.25 1]);
        set(gca,'YTickLabel','');
        set(gca,'FontSize',24);
        set(gca,'xlim',[0 size(CurMat,2)+1])
        title(int2nt(k))

        for x = 1:length(TextXpos)
            for y = 1:size(TextYpos,1)
                if CurMat(y,x) == 0
                    continue
                end
                Ytext = num2str(CurMat(y,x),'%0.2f');
                text(TextXpos(x), TextYpos(y,x), Ytext, 'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',12);
            end
        end
    end
end