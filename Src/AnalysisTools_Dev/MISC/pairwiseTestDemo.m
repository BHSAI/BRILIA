%pairwiseTestDemo demonstrates how pairwiseTest determines a p-value for
%the difference of two vectors of varying sizes, in comparison to the
%t-test. The t-test will reveal how different the 2 "distributions" are,
%but the don-test will reveal how well you can distinguish the 2 values as
%being different from a 0 median, based on all combination of differences.

function pairwiseTestDemo
X = randn(1, 50)+4;
Y = randn(1, 50)+6;
Alpha = 0.05;
Tail = 'L';
pairwiseTestDemoCalc(X, Y, Alpha, Tail)

function pairwiseTestDemoCalc(X, Y, Alpha, Tail)
[Ps, Hs] = signtest(X, Y, 'alpha', Alpha, 'Tail', Tail);
[Pr, Hr] = signrank(X, Y, 'alpha', Alpha, 'Tail', Tail);
[Ht, Pt, CIt] = ttest2(X, Y, 'alpha', Alpha, 'Tail', Tail);
[Hd, Pd, CId, DXY] = pairwiseTest(X, Y, 'alpha', Alpha, 'Tail', Tail);

E = [min([X(:); Y(:)])-0.5:0.3:max([X(:); Y(:)])+0.5];
Nx = histcounts(X, E);
Ny = histcounts(Y, E);

subplot(2, 1, 1)
Z = (E(1:end-1) + E(2:end)) / 2; 
Bx = bar(Z, [Nx(:), Ny(:)]);
Bx(1).FaceColor = 'k';
Bx(2).FaceColor = 'r';
legend('X', 'Y', 'Location','NorthEast')
set(gca, 'XTick', Z, 'XTickLabelRotation', 90)
title(['Want to test if X < Y at ' sprintf('Alpha = %0.2f', Alpha)]);
xlabel('Value')
ylabel('Freq.');

subplot(2, 1, 2)
[Nxy, Exy] = histcounts(DXY);
Zxy = (Exy(1:end-1)+Exy(2:end))/2;
Bx = bar(Zxy, Nxy);
Bx.FaceColor = 'k';
set(gca, 'XTick', Zxy, 'XTickLabelRotation', 90, 'XLim', [Exy(1), Exy(end)])
hold(gca, 'on')
if Ps/0.0001 < 1
    T{1} = sprintf('sign-test: H = %0.0f, p << 0.001', Hs);
else
    T{1} = sprintf('sign-test: H = %0.0f, p = %0.4f', Hs, Ps);
end

if Pr/0.0001 < 1
    T{2} = sprintf('signrank-test: H = %0.0f, p << 0.001', Hr);
else
    T{2} = sprintf('signrank-test: H = %0.0f, p = %0.4f', Hr, Pr);
end

if Pt/0.0001 < 1
    T{3} = sprintf('t-test: H = %0.0f, p << 0.001, CI = [%0.2f, %0.2f]', Ht, CIt(1), CIt(2));
else
    T{3} = sprintf('t-test: H = %0.0f, p = %0.4f, CI = [%0.2f, %0.2f]', Ht, Pt, CIt(1), CIt(2));
end

if Pd/0.0001 < 1
    T{4} = sprintf('don-test: H = %0.0f, p << 0.001, CI = [%0.2f, %0.2f]', Hd, CId(1), CId(2));
else
    T{4} = sprintf('don-test: H = %0.0f, p = %0.4f, CI = [%0.2f, %0.2f]', Hd, Pd, CId(1), CId(2));
end

title(T);
xlabel('Y_i - X_j');
ylabel('Freq.');
drawnow
plot([0 0], get(gca, 'ylim'), 'r');
hold(gca, 'off')
drawnow
resizeSubplots;