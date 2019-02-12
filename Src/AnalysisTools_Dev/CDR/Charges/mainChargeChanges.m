function mainChargeChanges(S)

for f = 1:length(S)
    S(f).ChargeChanges = calcChargeChanges(S(f).VDJdata, S(f).VDJheader);
end

for f = 1:length(S)
    Freq = vertcat(S(f).ChargeChanges.hCDR3{:, 2});
    N = histc(Freq, -5:5);
    
    figure(f)
    bar(-5:5, N);
    xlabel('Charge Change');
    ylabel('Freq. of Parent-Child CDR3 charge changes')
    TitleName = getTitleFileName(S(f).FileName);
    title(TitleName);
    resizeSubplots;
    %savePlot(gcf, 'saveas', [num2str(f) '_ChargeChangeDir.png']);
end