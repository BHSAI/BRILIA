function plotChargeChanges(S)

Freq = S.ChargeChanges.hCDR3;
Freq(Freq(:, 1) == 0, :) = [];
UnqLen = unique(Freq(:, 1));
for y = 1:length(UnqLen)
    Idx = UnqLen(y) == Freq(:, 1);
    N = histc(Freq(Idx, 2), [-10:10]);
    bar([-10:10], N, 'histc');
    title(num2str(UnqLen(y)));
end