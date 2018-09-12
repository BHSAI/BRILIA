%loadDemoPlot will creat a demo plot used to show how FigureTools work.
%
%  Gx = loadDemoPlot;
%
function Gx = loadDemoPlot()
X = 1:3:10;
Y1 = X.^2;
Y2 = X / 2;
Y3 = -X;
Gx = figure();
subplot(2,2,1) 
plot(X, Y1, 'r');
subplot(2,2,2) 
plot(X, Y2, 'g');
subplot(2,2,3) 
plot(X, Y3, 'b');
