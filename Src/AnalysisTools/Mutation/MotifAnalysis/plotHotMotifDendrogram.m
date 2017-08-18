%plotHotMotifDendrogram will plot dendrograms to compare similarity in
%mutation frequency distribution amongst hotspot motifs.
%
%  ImageName = plotHotMotifDendrogram(HotMotifData, Param, Value, ...)
%
%  INPUT
%    HotMotifData: output of getMotifData, which contains information on
%      hotspot motifs and pairwise mutations frequencies.
%    Param: parameter names of axes, such as 'XLim', 'YLim', 'Visible'.
%    Value: values for the parameters.
%
%  OUTPUT
%    ImageName: the full file name where the plot is saved.
%
%  NOTE
%    This function invokes resizeFigure, resizeSubplots and savePlots,
%    hence Param-Value pairs for these functions can be used here.

function varargout = plotHotMotifDendrogram(HotMotifData, varargin)
P = inputParser;
addParameter(P, 'FigWidth', 6, @(x) isnumeric(x) && x > 1);
addParameter(P, 'FigHeight', 3.5, @(x) isnumeric(x) && x > 1);
addParameter(P, 'FontName', 'Arial', @(x) ischar(x));
addParameter(P, 'FontSize', 12, @(x) isnumeric(x) && x > 1);
addParameter(P, 'SaveSubDir', 'Analysis', @ischar);
addParameter(P, 'DefaultSaveAs', 'HotMotifDendrogram.png', @ishcar);
addParameter(P, 'Visible', 'on', @(x) ischar(x) && ismember(lower(x), {'on', 'off'}));
[Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInput(P, varargin{:});
if ReturnThis
   varargout = {Ps, Pu, ExpPs, ExpPu};
   return;
end
FigWidth = Ps.FigWidth;
FigHeight = Ps.FigHeight;
FontName = Ps.FontName;
FontSize = Ps.FontSize;
SaveSubDir = Ps.SaveSubDir;
DefaultSaveAs = Ps.DefaultSaveAs;
Visible = Ps.Visible;

Gx = figure('Visible', Visible);
pause(0.05); %Wait 50ms for concurrency issues
resizeFigure(Gx, FigWidth, FigHeight, Visible);

Axs = gobjects(4, 1);
MaxXLim = [0 0];
HotIdxAdj = [-1 0;  %A is the 2nd letter, so it's 2+[-1]: 2+[0]
             -2 0;  %C is the 3rd letter, so it's 3+[-2]: 3+[0]
              0 2;  %G is the 1st letter, so it's 1+ [0]: 1+[2]
              0 1]; %T is the 1st letter, so it's 1+ [0]: 1+[1]
for f = 1:4
    %Get the data and normalize
    Data = cell2mat(HotMotifData.(int2nt(f))(:, 2:end));
    NormData = Data ./ repmat(sum(Data, 2), 1, 4); %Remember to transpose due to calcDistribDist requirements
    MotifList = HotMotifData.(int2nt(f))(:, 1);

    %Delete out motifs with 0 occurrences
    DelLoc = isnan(NormData(:, 1));
    MotifList(DelLoc) = [];
    NormData(DelLoc, :) = [];
    
    %Compute distribution distances between motifs
    DistMat = calcDistribDist(NormData', 'Method', 'bhatt');
    InfLoc = isinf(DistMat);
    if max(InfLoc(:)) > 0
        DistMat(InfLoc(:)) = 10*max(DistMat(~InfLoc(:)));
    end
    Link = linkage(DistMat, 'single');
    LeafOrder = optimalleaforder(Link, DistMat);

    if f == 1
        subplot(2, 3, 2)
    elseif f == 2
        subplot(2, 3, 1)
    elseif f == 3
        subplot(2, 3, 4)
    elseif f == 4
        subplot(2, 3, 3)
    end
    %Plot dendrogram onto separate plots
    Dx = dendrogram(Link, 'Reorder', LeafOrder, 'Orientation', 'left');
    %Dx = dendrogram(Link, 'Orientation', 'left');
    set(Dx, 'LineWidth', 1.5, 'Color', [0 0 0])
    XLim = get(gca, 'XLim');
    if XLim(end) > MaxXLim(end)
        MaxXLim(end) = XLim(end);
    end
    MotifLabel = MotifList(LeafOrder);
    for j = 1:length(MotifLabel)
        MutLoc = eval(MotifLabel{j}(4));
        MotifRange = MutLoc + HotIdxAdj(f, :);
        TriNuc = lower(MotifLabel{j}(1:3));
        TriNuc(MotifRange(1):MotifRange(end)) = upper(TriNuc(MotifRange(1):MotifRange(end)));
        MotifLabel{j} = TriNuc;
    end
    set(gca, 'YTickLabel', MotifLabel);
    Axs(f) = gca;
end

%Determine the standard X axis properties
MaxXLim(end) = ceil(MaxXLim(end) / 0.05) * 0.05;
XIncr = floor((MaxXLim(end) / 3) / 0.05) * 0.05;
XTick = [MaxXLim(1):XIncr:MaxXLim(end)];

%Format all plots to final form, except Y labels
HotSpotMotif = {'WA', 'WRC', 'GYW', 'TW'};
UnderlinePos = [2 3 1 1]; %Which letter of motif to underline
for f = 1:4
    set(Axs(f), 'FontName', 'Courier New', 'FontSize', FontSize, 'LineWidth', 1.5, 'Box', 'on');
    set(Axs(f), 'XTick', XTick, 'XLim', MaxXLim)
    title(Axs(f), getLatexUnderline(HotSpotMotif{f}, UnderlinePos(f)), 'Interpreter', 'Latex', 'FontName', FontName, 'FontSize', FontSize);
    setPlotTickDecimal(Axs(f), 2, -1); 
end

%Replace YTickLabel with text in order to get the formatted text option.
for f = 1:4
    YTickLabel = get(Axs(f), 'YTickLabel');
    for j = 1:length(YTickLabel)
        UpperCaseLoc = regexp(YTickLabel{j}, '[A-Z]');
        UnderlineLoc = UpperCaseLoc(UnderlinePos(f));
        UnderYTickLabel = getLatexUnderline(YTickLabel{j}, UnderlineLoc);
        YTickLabel{j} = UnderYTickLabel;% get(Axs(f), 'YTickLabel');
    end
    set(Axs(f), 'YTickLabel', YTickLabel, 'TickLabelInterpreter', 'Latex');
end

resizeSubplots(Gx, 'ScaleVertical', 'y', 'HorzSpacer', 0.01, 'VertSpacer', 0.01, 'FigSpacer', 0.01, ExpPu{:});

FullSaveName = prepSaveTarget(ExpPu{:}, 'SaveSubDir', SaveSubDir, 'SaveExt', '.png', 'DefaultSaveAs', DefaultSaveAs, 'MakeSaveDir', 'y');
savePlot(Gx, ExpPu{:}, 'SaveAs', FullSaveName);

if strcmpi(Visible, 'off')
    close(Gx);
end

varargout{1} = FullSaveName;
