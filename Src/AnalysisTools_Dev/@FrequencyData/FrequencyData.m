%FrequencyData can be plotted many ways, such as PDF, CDF, Confetti. It can
%also be used to plot diversity indices. Use this object to store data of
%frequencies of any entity.
%
%  Obj = FrequencyData(Data1, Data2, Data2, ..., Param, Value, ...);
%
%  INPUT
%    DataN: Px2 array where N(:,1) are values and N(:,2) are frequencies
%           OR Qx1 array where each element is a unique observation
%    Param-Value are for the setable fields of the object. Type 
%           properties('FrequencyData') to see a list of setable fields
%
%  METHODS
%    plotCDF  &  saveCDF
%    plotPDF  &  savePDF
%    plotConfetti   &  saveConfetti
%    plotDiversity  &  saveDiversity
%    plotEntropy    &  saveEntropy
%
%  EXAMPLE
%    Data1 = {}
%
classdef FrequencyData < DataInterface
    properties (Access = public)
        XTickLabel cell %Cell array for X-axis
    end
     
    properties (GetAccess = public, SetAccess = private)
        MeanHandle     matlab.graphics.Graphics %Mean handles
        ErrorHandle    matlab.graphics.Graphics %Mrror handles
        SigTestHandle  matlab.graphics.Graphics %Text handles of stat. sig. pairwise comparisons
        GroupStat      struct %Results of groupwiseTest
        Alpha          double %Alpha value for doing groupwise comparisons
        ErrorMode char = 'minmax' %Determine the error mode to use
    end
       
    %Methods for formatting and checking the input data
    methods (Static)
        function Data = format(varargin)
            %Formats input data such that it can be accepted by
            %DataInterface objects for this particular subclass.
            Data = cell(nargin, 1);
            for j = 1:nargin
                Data{j} = formatData(varargin{j});
            end
            
            function Data = formatData(Data)
                if min(size(Data)) == 1
                    Data = Data(:);
                end
                switch size(Data, 2)
                    case 1 %Needs to be counted into a Mx2 cell array
                        Data = countData(Data);
                    case 2
                        if ~iscell(Data)
                            Data = num2cell(Data);
                        elseif ~all(cellfun(@isdouble, Data(:, 2)))
                            error('%s: Input data 2nd column must be all double for frequencies.', mfilename);
                        end
                    otherwise
                        error('%s: Input data must be a vector of obersvations, or a Mx2 array of [ObservedValue Frequency]', mfilename);
                end
            end
        end
    end
    
    methods (Access = private)
        function Data = getDataAsCell(O, Option)
            N = numel(O.Data);
            Data = cell(1, N);
            if nargin < 2
                Option = '';
            end
            for j = 1:numel(O.Data)
                Freq = O.Data(j).Data;
                if all(cellfun(@isnumeric, Freq(:, 1))) %These are numbers that can be multiplied
                    NumVal = cell2mat(Freq(:, 1));
                    NumFreq = cell2mat(Freq(:, 2));
                    ExpFreq = repelem(NumVal, NumFreq);
                    Data{j} = sort(ExpFreq);
                else %Cannot plot categorical data
                    disp('Cannot plot CDF for categorical data frequencies.');
                    return
                end
            end
            
            DoPercentage = endsWith(Option, 'perc', 'ignorecase', true);
            DoNormalize = startsWith(Option, {'norm', 'normalize'}, 'ignorecase', true);
            for f = 1:numel(O.ModData)
                if DoNormalize
                    Data{f} = Data{f} ./ sum(Data{f});
                    if DoPercentage
                        Data{f} = Data{f} * 100;
                    end
                end
            end
        end
    end
    
    methods (Access = public)
        function O = FrequencyData(varargin)
            %Creates a FrequencyData object
            O = O@DataInterface; %Use superclass constructor (REQUIRED)
            O.setGx('PDF', []);
            O.setGx('CDF', []);
            O.setGx('Confetti', []);
            O.setGx('Diversity', []);
            
            %Extract and add the Data to object
            DataIdx = find(cellfun(@ischar, varargin), 1)-1;
            if isempty(DataIdx); DataIdx = numel(varargin); end
            O.Data = FrequencyData.format(varargin{1:DataIdx});
            varargin = varargin(DataIdx+1:end);
            
            %Parse the remaining inputs
            P = inputParser;
            P.addParameter('DataName', mfilename, @ischar);
            P.addParameter('Group', [], @isnumeric);
            P.addParameter('GroupName', {}, @isnumeric);
            P.addParameter('CtrlGroup',  1, @isnumeric);
            P.addParameter('Width',  6, @(x) isnumeric(x) && x > 0.1);
            P.addParameter('Height', 3, @(x) isnumeric(x) && x > 0.1);
            P.addParameter('Alpha', 0.05, @(x) isnumeric(x) && x > 0 && x < 1);
            P.parse(varargin{:});
            
            %Update the object based on parsed inputs
            O.DataName = P.Results.DataName;
            O.Group  = P.Results.Group;
            O.GroupName = P.Results.GroupName;
            O.CtrlGroup = P.Results.CtrlGroup;
            O.Width  = P.Results.Width;
            O.Height = P.Results.Height;
            O.Alpha  = P.Results.Alpha;
            
            O.reset;  %Resets the ModData field
        end
        
        function add(O, Data, varargin)
            %Adds a data point to FrequencyData
            Data = FrequencyData.format(Data);
            O.add@DataInterface(Data, varargin{:});
            O.reset;
        end
        
        function del(O, Number)
            %Deletes a data point from FrequencyData
            O.del@DataInterface(Number);
            O.reset;
        end
             
        function filter(O)
            %Filters the stored data according to some parameter-value pairs
            disp('There is no filter option for this Data.');
        end
        
        function pool(O)
            %Pools data by the assigned Group
            disp('"pool" has not been implemented.');
        end

        function reset(O)
            %Resets the FrequencyData's internal data structure
            
            %Extracts the frequency data from the raw data
            if numel(O.Data) == 1 && isempty(O.Data.Data)
                disp('There is no data left.');
                O.ModData = struct;
                return
            end
            ConformData = conformDist(O.Data.Data);
            O.ModData.Label = ConformData(:, 1);
            O.ModData.Indiv = cell2mat(ConformData(:, 2:end));
            O.ModData.MultFactor = ones(1, size(O.Data, 2));
            O.ModData.PairH = [];
            O.setError(O.ErrorMode);
        end
        
        function open(O)
            %Opens a csv file
            disp('"open" has not been implemented.');
        end
        
        function save(O)
            %Saves the data to a csv file
            if isempty(O.SaveDir)
                warning('%s: No save directory set.', mfilename);
                return
            end
            SaveHeader = ['Label', cellfun(@(x, y) sprintf('%d_Grp%d', x, y), num2cell(1:numel(O.Group)), num2cell(O.Group), 'un', 0)];
            SaveData = [O.ModData.Label num2cell(O.ModData.Indiv)];
            SaveName = fullfile(O.SaveDir, [O.DataName '.csv']);
            writeDlmFile([SaveHeader; SaveData], SaveName);
            fprintf('%s: Save to file "%s".\n', mfilename, SaveName);
        end
    end
    
    %Subclass special data-modifying methods
    methods
        function rescale(O)
            %Rescales the frequencies to match the mean of the CtrlGroup values
            Idx = find(O.Group == O.CtrlGroup);
            Idx = ternary(isempty(Idx), 1, Idx);
            [O.ModData.Indiv, O.ModData.MultFactor] = rescaleDist(O.ModData.Indiv, O.ModData.Indiv(:, Idx));
            O.ModData.PairH = [];
            O.setError;
        end

        function rebin(O, NewBin)
            %Rebins the data accordin to new specified bin values, keeping the mid point values
            if nargin < 2
                NewBin = [];
            end
            Data = num2cell(O.ModData.Indiv);
            Label = O.ModData.Label;
            NewData = rebinData([Label, Data], NewBin);
            O.ModData.Indiv = cell2mat(NewData(:, 2:end));
            O.ModData.Label = NewData(:, 1);
            O.ModData.PairH = [];
            O.setError;
        end
        
        function normalize(O)
            %Normalizes frequencies so that sum of frequencies = 1
            O.ModData.Indiv = O.ModData.Indiv./sum(O.ModData.Indiv);
            O.ModData.PairH = [];
            O.ModData.MultFactor = ones(1, size(O.Data, 2));
            O.setError;
        end
        
        function setError(O, Option)
            %Sets the error bars according to "minmax" of data or "std" around mean of data.
            if nargin == 2
                PrevOption = O.ErrorMode;
                O.ErrorMode = Option;
            end
            switch lower(O.ErrorMode)
                case 'minmax'
                    MaxData = applyall('max', O.ModData.Indiv, O.Group);
                    MinData = applyall('min', O.ModData.Indiv, O.Group);
                    MidData = (MaxData + MinData) / 2;
                    O.ModData.Error = MaxData - MidData;
                    O.ModData.Group = MidData;
                case 'std'
                    O.ModData.Group = applyall('mean', O.ModData.Indiv, O.Group);
                    O.ModData.Error = applyall('std',  O.ModData.Indiv, O.Group);
                otherwise
                    warning('%s: Invalid Error Mode option. Use ''minmax'' or ''std''. Using previous mode.', mfilename);
                    O.ErrorMoode = PrevOption;
            end
        end
        
        function groupwiseTest(O, varargin)
            %Conducts a group-wise comparison of frequencies to determine siginificant difference per bar (PDF)
            [~, Gall]= groupwiseTest(O.ModData.Indiv, O.Group, 'Alpha', O.Alpha, 'Size', cellfun(@(x) sum(cell2mat(x(:, 2))), {O.Data.Data}), varargin{:});
            O.ModData.PairH = Gall(:, 3:end)' <= O.Alpha;
        end
    end
    
    %Plot methods
    methods
        function plot(O)
            disp('This has not been coded yet');
        end
        
        function plotPDF(O, TopN)
            if ~O.hasData
                disp('There is no data left to plot.');
                return
            end
            if ~O.isGxValid('PDF')
                O.Gx.PDF = figure;
            end
            cla(O.Gx.PDF);

            %Plots all data or the top N bars
            if nargin > 1 && size(O.ModData.Indiv, 1) >= TopN
                DataSum = max(O.ModData.Indiv, [], 2);
                [~, Idx] = sort(DataSum, 'descend');
                ShowIdx = sort(Idx(1:TopN));
            else
                ShowIdx = [];
            end
            %Plot the bar graph
            [~, O.MeanHandle, O.ErrorHandle, O.SigTestHandle] = plotErrorBox(O.Gx.PDF, O.ModData.Group, 'STD', O.ModData.Error, 'PairH', O.ModData.PairH, 'XTickLabel', O.ModData.Label, 'IndivData', O.ModData.Indiv, 'Group', O.Group, 'ShowIdx', ShowIdx);
            if ~isempty(O.DataName)
                TitleName = strrep(sprintf('%s', strrep(O.DataName, '_', ' ')), '_', '\_');
                setAxes(O.Gx.PDF, 'Title-String', TitleName);
            end
            O.applyDefault('PDF');
        end
        
        function plotCDF(O)
            if ~O.hasData
                disp('There is no data left to plot.');
                return
            end
            Data = O.getDataAsCell;
            if all(cellfun('isempty', Data))
                return
            end
            if ~O.isGxValid('CDF')
                O.Gx.CDF = figure;
            end
            N = numel(O.Data);
            StdColor = getStdColor;
            hold(O.Ax.CDF, 'on')
            MaxX = 0;
            for j = 1:N
                X = Data{j}(:);
                Y = [1:numel(X)]/numel(X);
                plot(O.Ax.CDF, X, Y, 'Color', StdColor(O.Data(j).Group, :), 'LineWidth', 1.5);
                if max(X) > MaxX
                    MaxX = max(X);
                end
            end
            hold(O.Ax.CDF, 'on')
            xlim(O.Ax.CDF, [0, MaxX])
            ylim(O.Ax.CDF, [0, 1]);
            setPlotTickDecimal(O.Ax.CDF, 0, 1)
            setAxes(O.Ax.CDF, 'XLabel-String', 'Value', 'YLabel-String', 'CDF')
            O.applyDefault('CDF')
        end
%             
%             
%             % LxCDF = gobjects(N, 1);
%             % LxYIntercept = gobjects(N, 1);
%             CDFs = zeros(N, 1);
%             
%             
%             Cutoff = 5;
%             
%             for j = 1:N
%                 Template = sort(O.ModData(j).Data);
%                 CDF = (1:numel(Template))/numel(Template);
%                 plot(O.Ax.CDF, Template, CDF, 'Color', StdColor(O.ModData(j).Group, :), 'LineWidth', 1.5);
%                 if j == 1
%                     hold(gca, 'on')
%                 end
%                 
%                 %Computer the y-intercept line
%                 YIdx = find(Template <= Cutoff, 1, 'last');
%                 CDFs(j) = CDF(YIdx);
%                 %    LxYIntercept(j) = scatter(O.Ax.CDF, Template(YIdx), CDFs(j) , 'o', 'MarkerEdgeColor', StdColor(O.Group(j), :));
%             end
%             % plot([Cutoff Cutoff], [0 1], 'k--');
%             hold(gca, 'off')
%             xlim([0 50])
%             
%             Group = [O.ModData.Group]';
%             
%             Mean = splitapply(@mean, CDFs, Group);
%             Std = splitapply(@std, CDFs, Group);
%             Entropy = zeros(N, 1);
%             for j = 1:numel(Entropy)
%                 [~, Entropy(j)] = calcDiversity('Shannon', O.ModData(j).Data);
%             end
%             AvgEntropy = splitapply(@mean, Entropy, Group);
%             StdEntropy = splitapply(@std, Entropy, Group);
%             
%             TotalS = cellfun('length', {O.ModData.Data})';
%             Evenness = Entropy ./ log(TotalS);
%             AvgEvenness = splitapply(@mean, Evenness, Group);
%             StdEvenness = splitapply(@std, Evenness, Group);
%             
%             %Assemble text legend
%             DataHeader = {'Group' sprintf('Fr. %s Cutoff (%d)', char(8804), Cutoff) 'Shannon Entropy', 'Evenness'};
%             
%             %Final Adjustements
%             setPlotTickDecimal(O.Ax.CDF, 0, 1);
%             title(O.Ax.CDF, [O.DataName 'CDF']);
%             resizeSubplots(O.Ax.CDF);
%             return
%             
%             GroupTxt = cellfun(@(x) sprintf('%d', x), num2cell(unique(O.Group)), 'un', 0);
%             CutoffTxt = cellfun(@(x,y) sprintf('%0.3f %s %0.3f', x, char(177), y), num2cell(Mean), num2cell(Std), 'un', 0);
%             ShannonTxt = cellfun(@(x,y) sprintf('%0.1f %s %0.1f', x, char(177), y), num2cell(AvgEntropy), num2cell(StdEntropy), 'un', 0);
%             EvennessTxt = cellfun(@(x,y) sprintf('%0.3f %s %0.3f', x, char(177), y), num2cell(AvgEvenness), num2cell(StdEvenness), 'un', 0);
%             DataValue = [GroupTxt, CutoffTxt, ShannonTxt, EvennessTxt];
%             DataTxt = [DataHeader; DataValue];
%             
%             MaxCharPerCol = max(cellfun('length', DataTxt), [], 1);
%             MaxCharPerCol = repmat(MaxCharPerCol, size(DataTxt, 1), 1);
%             for j = 1:numel(DataTxt)
%                 DataTxt{j} = sprintf('%s%s', DataTxt{j}, repelem(' ', 1, MaxCharPerCol(j) - numel(DataTxt{j})));
%             end
%             
%             LegendTxt = cell(size(DataTxt, 1), 1);
%             for j = 1:numel(LegendTxt)
%                 LegendTxt{j} = sprintf('%s  ', DataTxt{j, :});
%             end
%             
%             X = Cutoff + 1;
%             Y = 0.5;
%             text(O.Ax.CDF, X, Y, LegendTxt, 'FontName', 'Courier New', 'FontSize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
%             
%             
%             %Final Adjustements
%             setPlotTickDecimal(O.Ax.CDF, 0, 1);
%             title(O.Ax.CDF, [O.DataName 'CDF']);
%             resizeSubplots(O.Ax.CDF);       
        function plotConfetti(O, varargin)
            Data = O.getDataAsCell;
            if all(cellfun('isempty', Data))
                return
            end
            P = inputParser;
            P.addParameter('DownSample', 'n', @(x) ismember(lower(x), {'y', 'n'}));
            P.parse(varargin{:});
            DownSample = strcmpi(P.Results.DownSample, 'y');

            %Down sample to the minimum count of all Data
            if DownSample
                DataCt = cellfun('length', Data);
                MinCt = min(DataCt);
                for j = 1:numel(Data)
                    Data{j} = sort(Data{j});
                    GetIdx = round(linspace(1, DataCt(j), MinCt));
                    Data{j} = Data{j}(GetIdx);
                end
            end
            
            %Make sure the confetti map matrices are made
            if ~isfield(O.ModData, 'Confetti') || isempty(O.ModData.Confetti)
                O.ModData.Confetti = cell(1, numel(Data));
            end
            %Calc and plot Confettis
            for j = 1:numel(Data)
                if isempty(O.ModData.Confetti{j}) || ...
                   size(O.ModData.Confetti{j}, 1) ~= size(O.ModData.Confetti{j}, 2) || ...
                   max(O.ModData.Confetti{j}(:)) ~= numel(Data{j})
                    fprintf('%s: Making ConfettiMap for Data # %d.\n', mfilename, j)
                    O.ModData.Confetti{j} = makeConfettiMap(Data{j});
                end
                FigName = sprintf('Confetti_%d', j);
                if ~O.isGxValid(FigName)
                    O.setGx(FigName, figure);
                end
                plotConfetti(O.Ax.(FigName), O.ModData.Confetti{j})
            end
            %Do not apply defaults for Confetti plots since you want the
            %WxH to be a square.
        end
        
        function plotEntropy(O, varargin)
            Data = O.getDataAsCell;
            if all(cellfun('isempty', Data))
                return
            end
            if isempty(varargin)
                varargin = {'richness'};
            end
            Entropy = zeros(numel(varargin), numel(Data));
            DelThis = zeros(numel(varargin), 1, 'logical');
            for j = 1:numel(varargin)
                try
                    for k = 1:numel(Data)
                        [~, Entropy(j, k)] = calcDiversity(varargin{j}, Data{k});
                    end
                catch
                    warning('%s: Unrecognized diversity option, "%s".', mfilename, varargin{1});
                    DelThis(j) = 1;
                    continue
                end
            end
            Entropy = Entropy(~DelThis, :);
            varargin = varargin(~DelThis);
            
            Mean = applyall('mean', Entropy, O.Group);
            Error = applyall('std', Entropy, O.Group);
            Ax = O.Ax.Simpson;
            if isempty(Ax) || ~isvalid(Ax)
                Ax = [];
            end
            Gx = plotErrorBox(Ax, Mean, 'Std', Error, 'IndivData', Entropy, 'Group', O.Group);
            Ax = get(Gx, 'children');
            O.Ax.Entropy = Ax;
            ylabel(Ax, 'Entropy');
            Ax.XTickLabel = varargin;
            resizeSubplots(Ax);            
        end
        
        function plotDiversity(O, varargin)
            Data = O.getDataAsCell;
            if all(cellfun('isempty', Data))
                return
            end
            if isempty(varargin)
                varargin = {'richness'};
            end
            Diversity = zeros(numel(varargin), numel(Data));
            DelThis = zeros(numel(varargin), 1, 'logical');
            for j = 1:numel(varargin)
                try
                    for k = 1:numel(Data)
                        Diversity(j, k) = calcDiversity(varargin{j}, Data{k});
                    end
                catch
                    warning('%s: Unrecognized diversity option, "%s".', mfilename, varargin{1});
                    DelThis(j) = 1;
                    continue
                end
            end
            Diversity = Diversity(~DelThis, :);
            varargin = varargin(~DelThis);
            
            Mean = applyall('mean', Diversity, [O.Group]);
            Error = applyall('std', Diversity, [O.Group]);
            if ~O.isAxValid('Diversity')
                O.setAx('Diversity', newplot);
            end
            Ax = O.Ax.Diversity;
            if isempty(Ax) || ~isvalid(Ax)
                Ax = [];
            end
            Gx = plotErrorBox(Ax, Mean, 'Std', Error, 'IndivData', Diversity, 'Group', O.Group);
            Ax = get(Gx, 'children');
            O.Ax.Diversity = Ax;
            ylabel(Ax, 'Diversity');
            Ax.XTickLabel = varargin;
            resizeSubplots(Ax);            
        end
        
        function plotTemplate(O, varargin)
            P = inputParser;
            addParameter(P, 'Option', 'normperc', @(x) ismember(lower(x), {'', 'norm', 'normalize', 'normperc'}));
            addParameter(P, 'YScale', 'linear', @(x) ismember(lower(x), {'log', 'linear'}));
            parse(P, varargin{:})
            Option = P.Results.Option;
            YScale = P.Results.YScale;
            
            if isempty(O.Ax.Template) || ~isvalid(O.Ax.Template)
                Gx = figure;
                O.Ax.Template = axes(Gx);
            end
            Ax = O.Ax.Template;
            
            Data = O.getDataAsCell(Option);
            for f = 1:numel(Data)
                N = numel(Data{f});
                scatter(Ax, f+0.8*rand(N, 1)-0.5, Data{f}, 10, O.Color(O.Group(f), :));
                if f == 1
                    hold(Ax, 'on')
                end
            end
            hold(Ax, 'off')
            Ax.YScale = YScale;
            if strcmpi(Option, 'norm')
                ylabel(Ax, 'Clonal Freq. Fr.')
                setPlotTickDecimal(Ax, -1, 3);
            elseif strcmpi(Option, 'normperc')
                ylabel(Ax, 'Clonal Freq. (%)')
                setPlotTickDecimal(Ax, -1, 1);
            else
                ylabel(Ax, 'Clonal Freq.')
                setPlotTickDecimal(Ax, -1, 0);
            end
            setAxes(Ax, 'XTick', 1:numel(Data), 'XTickLabel', O.Group, 'Title-String', O.DataName);
            resizeSubplots(Ax)
        end
        
        %Template distribution where each bar is a unique entity, and each Y is the normalized probability
        %This shows the evenness of distributions well.
        function plotExpandedTemplate(O, varargin)
            P = inputParser;
            addParameter(P, 'Option', 'normperc', @(x) ismember(lower(x), {'', 'norm', 'normalize', 'normperc'}));
            addParameter(P, 'YScale', 'linear', @(x) ismember(lower(x), {'log', 'linear'}));
            parse(P, varargin{:})
            Option = P.Results.Option;
            YScale = P.Results.YScale;
            
            if ~isfield(O.Ax, 'ExpTempmlate') || isempty(O.Ax.ExpTemplate) || ~isvalid(O.Ax.ExpTemplate)
                Gx = figure;
                O.Ax.ExpTemplate = axes(Gx);
            end
            Ax = O.Ax.ExpTemplate;

            Data = O.getDataAsCell(Option);
            for f = 1:numel(Data)
                Data{f} = sort(Data{f}, 'descend');
                X = [1:numel(Data{f})]/numel(Data{f});
                if strcmpi(Option, 'normperc')
                    Y = Data{f}/sum(Data{f}) * 100; %in percentage
                else
                    Y = Data{f}/sum(Data{f}); %in percentage
                end
                plot(Ax, X, Y, 'Color', O.Color(O.Group(f), :));
                if f == 1
                    hold(Ax, 'on')
                end
            end
            hold(Ax, 'off')
            Ax.YRuler.Exponent = 0;
            setPlotTickDecimal(Ax, 1, 2)
            Ax.YScale = YScale;
            if strcmpi(Option, 'norm')
                ylabel(Ax, 'Clonal Freq. Fr.')
                setPlotTickDecimal(Ax, -1, 3);
            elseif strcmpi(Option, 'normperc')
                ylabel(Ax, 'Clonal Freq. (%)')
                setPlotTickDecimal(Ax, -1, 1);
            else
                ylabel(Ax, 'Clonal Freq.')
                setPlotTickDecimal(Ax, -1, 0);
            end
            setAxes(Ax, 'XTick', 1:numel(Data), 'XTickLabel', O.Group, 'Title-String', O.DataName);
            resizeSubplots(Ax)
        end
        
        
        
        
        function savePlot(O)
            %Saves all plot to the path "SaveDir\DataName.png"
            if isempty(O.SaveDir)
                warning('%s: No save directory set.', mfilename);
                return
            end
            if ~isempty(O.Gx) && isvalid(O.Gx)
                SaveName = fullfile(O.SaveDir, [O.DataName '.png']);
                savePlot(O.Gx, 'SaveAs', SaveName, 'DPI', 600);
            else
                warning('%s: No figure handle found.', mfilename);
                return
            end
        end
        
        function colorBars(O, RGB)
            %Colors each main bar a certain color
            arrayfun(@(x) set(x, 'FaceColor', RGB), O.MeanHandle);
        end
 
        
        function showErrorBars(O)
            %Shows error bars
            arrayfun(@(x) set(x, 'Visible', 'on'), O.ErrorHandle);
        end

        function hideErrorBars(O)
            %Hides error bars
            arrayfun(@(x) set(x, 'Visible', 'off'), O.ErrorHandle);
        end

        function colorErrorBars(O, RGB)
            %Colors each error bar a certain color
            arrayfun(@(x) set(x, 'FaceColor', RGB), O.ErrorHandle);
        end
        
        function showBars(O)
            %Shows main bars
            arrayfun(@(x) set(x, 'Visible', 'on'), O.MeanHandle);
        end
        
        function hideBars(O)
            %Hides main bars
            arrayfun(@(x) set(x, 'Visible', 'off'), O.MeanHandle);
        end
        
        function colorSigTest(O, RGB)
            %Colors each error bar a certain color
            arrayfun(@(x) set(x, 'Color', RGB), O.SigTestHandle);
        end
        
        function showSigText(O)
            %Shows main bars
            arrayfun(@(x) set(x, 'Visible', 'on'), O.SigTestHandle);
        end
        
        function hideSigText(O)
            %Hides main bars
            arrayfun(@(x) set(x, 'Visible', 'off'), O.SigTestHandle);
        end
    end
end