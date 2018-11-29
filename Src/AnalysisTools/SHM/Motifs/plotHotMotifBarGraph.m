%plotHotMotifBarGraph will plot the relative A, C, G, T mutations per each
%hotspot motif.
%
%  ImageName = plotHotMotifBarGraph(HotMotifData, Param, Value, ...)
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
%    This function invokes setAxes, resizeSubplots, and savePlots, hence
%    Param-Value pairs for these functions can be used here.
function varargout = plotHotMotifBarGraph(HotMotifData, varargin)
P = inputParser;
addParameter(P, 'FigWidth', 6, @(x) isnumeric(x) && x > 1);
addParameter(P, 'FigHeight', 6, @(x) isnumeric(x) && x > 1);
addParameter(P, 'YLim', [0 1.2], @(x) isnumeric(x) && length(x) == 2);
addParameter(P, 'XLim', [0.5 4.5], @(x) isnumeric(x) && length(x) == 2);
addParameter(P, 'XTickLabel', {'A', 'C', 'G', 'T'}, @(x) iscell(x) && length(x) == 4);
addParameter(P, 'FaceColor', [0 0 0], @(x) isnumeric(x) && length(x) == 3);
addParameter(P, 'FontName', 'Arial', @(x) ischar(x));
addParameter(P, 'FontSize', 12, @(x) isnumeric(x) && x > 1);
addParameter(P, 'SaveSubDir', 'Analysis', @ischar);
addParameter(P, 'DefaultSaveAs', 'HotMotifBarGraph.csv', @ishcar);
addParameter(P, 'Visible', 'on', @(x) ischar(x) && ismember(lower(x), {'on', 'off'}));
[Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInput(P, varargin{:});
if ReturnThis
   varargout = {Ps, Pu, ExpPs, ExpPu};
   return;
end
FigHeight = Ps.FigHeight;
FigWidth = Ps.FigWidth;
XLim = Ps.XLim;
YLim = Ps.YLim;
XTickLabel = Ps.XTickLabel;
FaceColor = Ps.FaceColor;
FontName = Ps.FontName;
FontSize = Ps.FontSize;
SaveSubDir = Ps.SaveSubDir;
DefaultSaveAs = Ps.DefaultSaveAs;
Visible = Ps.Visible;

HotIdxAdj = [-1 0;  %A is the 2nd letter, so it's 2+[-1]: 2+[0]
             -2 0;  %C is the 3rd letter, so it's 3+[-2]: 3+[0]
              0 2;  %G is the 1st letter, so it's 1+ [0]: 1+[2]
              0 1]; %T is the 1st letter, so it's 1+ [0]: 1+[1]
UnderlineAt = [2 3 1 1]; %Underlines hotspot at this position for 'W[A]<' 'WR[C]' '[G]YR' '[T]A'};

%There are 12 plots, so want to organized as C[4], G[4], A/T[4] mutations
Gx = figure('Visible', Visible);
drawnow('nocallbacks');
resizeFigure(Gx, 'FigWidth', FigWidth, 'FigHeight', FigHeight, 'Visible', Visible);

SubPlotNum = [1 1 1]; % for C, G, A/T
for f = 1:4
    %Determine unique hotspot motifs (A/T need reduction from 3 to 2 nts)
    Data = cell2mat(HotMotifData.(int2nt(f))(:, 2:end));
    TriNuc = HotMotifData.(int2nt(f))(:, 1);
    MutLoc = zeros(size(TriNuc));
    MotifNuc = cell(size(TriNuc));
    for j = 1:length(MutLoc)
        MutLoc = convStr2NumMEX(TriNuc{j}(4));
        MotifRange = MutLoc + HotIdxAdj(f, :);
        MotifNuc{j} = TriNuc{j}(MotifRange(1):MotifRange(end));
    end
    
    %Sum up the mut counts per each unique motif
    [UnqMotifNuc, ~, UnqIdx] = unique(MotifNuc);
    UnqMotifData = zeros(length(UnqMotifNuc), 4);
    for j = 1:max(UnqIdx)
        Idx = j == UnqIdx;
        UnqMotifData(j, :) = sum(Data(Idx, :), 1);
    end
    
    %Normalize the data
    UnqMotifData = UnqMotifData ./ repmat(sum(UnqMotifData, 2), 1, 4);
        
    %Plot into subplots
    for j = 1:length(UnqMotifNuc)
        if f == 1 || f == 4 %A/T
            s = 3 * SubPlotNum(3);
            SubPlotNum(3) = SubPlotNum(3) + 1;
        elseif f == 2 %C
            s = 3 * (SubPlotNum(1) - 1) + 1;
            SubPlotNum(1) = SubPlotNum(1) + 1;
        else
            s = 3 * (SubPlotNum(2) - 1) + 2;
            SubPlotNum(2) = SubPlotNum(2) + 1;
        end
        subplot(4, 3, s)
        Bx = bar(UnqMotifData(j, :));
        set(Bx, 'FaceColor', FaceColor);
        set(gca, 'FontName', FontName, 'FontSize', FontSize, 'XLim', XLim, 'YLim', YLim, 'XTickLabel', XTickLabel);
        setPlotTickDecimal(gca, -1, 1)
        setAxes(gca, ExpPu{:});
        TitleName = getLatexUnderline(UnqMotifNuc{j}, UnderlineAt(f));
        text(sum(XLim)/2, YLim(2) - 0.05, TitleName, 'VerticalAlignment', 'Top', 'HorizontalAlignment', 'Center', 'FontName', 'Arial', 'FontSize', 12, 'Interpreter', 'latex');
        if s < 10
            set(gca, 'XTickLabel', '');
        end
    end
end

resizeSubplots(Gx, 'ScaleVertical', 'y', 'HorzSpacer', 0.01, 'VertSpacer', 0.01, 'FigSpacer', 0.01, ExpPu{:});

FullSaveName = prepSaveTarget(ExpPu{:}, 'SaveSubDir', SaveSubDir, 'SaveExt', '.png', 'DefaultSaveAs', DefaultSaveAs, 'MakeSaveDir', 'y');
savePlot(Gx, ExpPu{:}, 'SaveAs', FullSaveName);
close(Gx);
varargout{1} = FullSaveName;
