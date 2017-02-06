function Gx = plotHotSpots(Amat,Cmat,Gmat,Tmat)

%To simplify the plotting, rearrange matrices
HotMat = {Amat{1} Cmat{1} Gmat{1} Tmat{1}};
AllMat = {Amat{2} Cmat{2} Gmat{2} Tmat{2}};
ColNum = size(Amat{1},2);

%Plotting all hot spots
for g = 1:2
    Gx{g} = figure;
    for k = 1:4
        subplot(1,4,k)
        if g == 1
            CurMat = HotMat{k};
            CurMat = CurMat / max(CurMat(:,3));
        else
            CurMat = AllMat{k};
            CurMat = CurMat / max(CurMat(:,3));
        end

        %Calculate the position for the texts
        TextXpos = [1:ColNum];
        TextYpos = cumsum(CurMat,1) - CurMat/2;

        H1 = bar(CurMat','stack');
        H1(1).BarWidth = 1;
        H1(1).BarWidth = 1;
        H1(1).FaceColor = [1.0 0.2 0.2];
        H1(2).FaceColor = [0.2 1.0 0.2];
        H1(3).FaceColor = [0.4 0.4 1.0];
        H1(4).FaceColor = [0.7 0.7 0.7];

        BarWidth = H1(1).BarWidth;
        ylim([0,1])
        xlim([TextXpos(1)-BarWidth/2 TextXpos(end)+BarWidth/2])
        Xlab = cell(1,ColNum);
        MidVal = ceil(ColNum/2);
        for q = 1:length(Xlab)
            Xlab{q} = sprintf('%d',q-MidVal);
        end        
        set(gca,'XTickLabel',Xlab)
        set(gca,'YTickLabel','');
        set(gca,'FontSize',24);
        set(gca,'OuterPosition',[0.25*(k-1) 0 0.25 1]);
        title(int2nt(k))

        for x = 1:length(TextXpos)
            for y = 1:size(TextYpos,1)
                if CurMat(y,x) == 0
                    continue
                end
                Ytext = sprintf('%2.0f',CurMat(y,x)*100); % num2str(CurMat(y,x),'%0.2f');
                text(TextXpos(x), TextYpos(y,x), Ytext, 'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',12,'FontName','Arial');
            end
        end
    end
end
