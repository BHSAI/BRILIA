%plotHotspotCount

function varargout = plotHotspotMutData(HotspotMutData)
%Compute the fractions of hotspot in cdr/total hotspots avail
Hcdr = HotspotMutData.Hcdr;
Hfwr = HotspotMutData.Hfwr;
Ncdr = HotspotMutData.Ncdr;
Nfwr = HotspotMutData.Nfwr;
ActualCDR = sum(Ncdr) / (sum(Ncdr) + sum(Nfwr));
ExpectedCDR = sum(Hcdr) / (sum(Hcdr) + sum(Hfwr));

%Group data analysis
GrpData = HotspotMutData.GrpData;
Y1 = GrpData(:, 2) ./ sum(GrpData(:, 2:3), 2);
Y2 = GrpData(:, 4) ./ sum(GrpData(:, 4:5), 2);
X = GrpData(:, 1);

Gx = figure();
subplot(2, 1, 1)
scatter(X, Y1, 10, [0.3 0.3 0.3])
Axs{1} = gca;
xlabel('Clonotype Size by Unique Seq')
ylabel('HotMut_{CDR} / HotMut_{All}')
title(sprintf('Actual / Exp. CDR3 Mut = %0.2f / %0.2f', ActualCDR, ExpectedCDR));
XLim = get(gca, 'xlim');
ylim([0 1]);
hold(gca, 'on')
scatter(X, Y2, 10, [1 0 0], 'filled')
legend('Expected', 'Actual')

subplot(2, 1, 2)
scatter(X, Y2 ./ Y1, 20, [0 0 0], 'filled')
ylim([0 20])
Axs{2} = gca;
xlabel('Clonotype Size by Unique Seq')
ylabel('Actual / Expected')
hold(gca, 'on')
plot([0 XLim(2)], [1 1], 'r')


if nargout >= 1
    varargout{1} =  Gx;
    if nargout >= 2
        varargout{2} = Axs;
    end
end
