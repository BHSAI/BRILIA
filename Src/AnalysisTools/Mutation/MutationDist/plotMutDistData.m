%plotMutDistData plots the mutation distance from the CDR3 ends, separated
%by the unique V gene family.

function [Gx, Ax] = plotMutDistData(MutDistData)
%CDR1 loc is the 27 to 38 AA residues
CDR1start = 26*3 + 1;
CDR1end = 38*3;

%CDR2 loc is the 56 to 65 AA residues
CDR2start = 55*3 + 1;
CDR2end = 65*3;

%CDR3 loc is the 105 AA residue
CDR3start = 104*3 + 1;

Vnames = MutDistData.V(:, 1);
Vdata = cell2mat(MutDistData.V(:, 2:end));
Vsum = sum(Vdata, 1);

Jnames = MutDistData.J(:, 1);
Jdata = cell2mat(MutDistData.J(:, 2:end));
Jsum = sum(Jdata, 1);

Gx = figure;
Ax = axes;

VJsum = [fliplr(Vsum) zeros(1, 15*3) Jsum];
VJtickLabel = [-length(Vsum):1 zeros(1, 15*3) 1:length(Jsum)];
VJtick = [1:length(VJsum)];
Bx = bar(VJsum);
Bx.FaceColor = [0 0 0];
xlim([0 length(VJsum)+1])
ylim([0 max(VJsum)]);
hold(gca, 'on');
plot([CDR1start CDR1end], max(VJsum)*0.95 * [1 1], 'r', 'LineWidth', 3);
plot([CDR2start CDR2end], max(VJsum)*0.95 * [1 1], 'b', 'LineWidth', 3);
plot([length(Vsum)+1 length(Vsum)+3*15-1], max(VJsum)*0.95 * [1 1], 'Color', [0 0.8 0], 'LineWidth', 3);

resizeFigure(Gx, 6, 3);
resizeSubplots(Gx, 'HorzSpacer', 0.05)
savePlot(Gx, 'Save', 'y', 'SaveAs', 'plotMutDistData.png')
% 
% 
% subplot(1, 2, 1)
% Bx1 = bar(fliplr(Vsum));
% title('V gene mutations from CDR3s')
% Bx1.FaceColor = [0 0 0];
% ylim([0 max(Vsum)]);
% hold(gca, 'on');
% plot([CDR1start CDR1end], max(Vsum)*0.95 * [1 1], 'r', 'LineWidth', 1);
% plot([CDR2start CDR2end], max(Vsum)*0.95 * [1 1], 'b', 'LineWidth', 1);
% 
% subplot(1, 2, 2)
% title('J gene mutations from CDR3s')
% Bx2 = bar(Jsum);
% Bx2.FaceColor = [0 0 0];
% ylim([0 max(Jsum)]);
% hold(gca, 'on');
% 
% resizeSubplots(Gx)
