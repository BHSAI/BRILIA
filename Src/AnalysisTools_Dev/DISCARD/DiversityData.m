%DiversityData is an object container for many different plots of templates. It
%handles the storage, manipulation, ploting, and saving of such data.
%
%  Methods
%    plot: plot ALL types of diversity plots for the current working data set
%    sort(varargin): sorts the x-axis label
%    viewTop: view top N bars, by mean value
%    viewBot: view top N bars, by mean value
%    savePlot: saves the current plot
%    saveData: saves the current data
%
%  Properties
%    Gx: cell array of figure handles
%    Mx: cell array of MxN box objects of means
%    Ex: cell array of MxN box objects of error bars
%    Hx: cell array of '*' objects for statistically different
%
%  NOTE: This is done per ALL REPERTOIRES
classdef DiversityData < DataInterface
%     properties (Dependent)
%         Group double    %grouping of data
%     end
    
    properties (Access = public)
        CtrlGroup %Ctrl group
        Level %Clustering level of template: Clone, Clonotype, CDR3, Kmer, etc
        ConfettiMap cell %A cell of MxM Confetti maps
        Color %Color map Mx3 matrix
%         Width          double %width in inch of figure
%         Height         double %height in inch of figure
    end
     
    properties (GetAccess = public, SetAccess = private)
        GroupStat      struct %stores information about the groupwise test
%         PrivateGroup   double %stores the current grouping of data, workaround dependency issues in MATLAB
    end
    
    methods (Access = public)
        function O = DiversityData(varargin)
            P = inputParser;
            addRequired(P, 'Data', @(x) iscell(x) || isnumeric(x));
            addParameter(P, 'DataName', mfilename, @ischar);
            addParameter(P, 'Group', [], @isnumeric);
            addParameter(P, 'CtrlGroup', [], @isnumeric);
            addParameter(P, 'Width',  6, @(x) isnumeric(x) && x > 0.1);
            addParameter(P, 'Height', 3, @(x) isnumeric(x) && x > 0.1);
            addParameter(P, 'Color', getStdColor, @(x) isnumeric(x) && size(X, 2) == 3);
            parse(P, varargin{:});
            O.Width = P.Results.Width;
            O.Height = P.Results.Height;
            O.Color = P.Results.Color;
            O.Data = P.Results.Data;
            O.DataName = P.Results.DataName;
            O.format; %Do this BEFORE assigning group, otherwise error can occur.
            O.Group = P.Results.Group;
            O.CtrlGroup = P.Results.CtrlGroup;
            O.reset;  %Resets the ModData field
            
            O.Ax = struct('Template', [], 'Diversity', [], 'Confetti', []);
        end
        
        function TF = check(O)
            TF = iscell(O.Data) && ~isempty(O.Data);
        end
        
        function filter(O, Option)
            if strcmpi(Option, 'pool') %pool data by group
                disp('Pooling data by group.');
                [UnqGroup, ~, UnqIdx] = unique(O.Group);
                O.ModData = struct; 
                for j = 1:numel(UnqGroup)
                    O.ModData(j).Data = vertcat(O.Data(UnqIdx == j).Data);
                    O.ModData(j).Idx = vertcat(O.Data(UnqIdx == j).Idx);
                    O.ModData(j).Group = UnqGroup(j);
                end
            elseif strcmpi(Option, 'downsample')
                MinSize = min(cellfun('length', {O.ModData.Data}));
                for j = 1:numel(O.ModData)
                    [Data, SortIdx] = sort(O.ModData(j).Data, 'descend');
                    GetIdx = round(linspace(1, numel(Data), MinSize));
                    O.ModData(j).Data = O.ModData(j).Data(SortIdx(GetIdx));
                    O.ModData(j).Idx = O.ModData(j).Idx(SortIdx(GetIdx));
                end
            end
        end
        
        function format(O)
            %All data should be of Mx1 cell of struct of .Data, .Idx. This will format everything into a single nonscalar struct.
            if isstruct(O.Data) && isfield(O.Data, 'Data') && isfield(O.Data, 'Idx')
                return
            elseif iscell(O.Data) && all(cellfun('isclass', O.Data, 'struct'))
                S(1:numel(O.Data)) = struct('Data', [], 'Idx', []);
                for k = 1:numel(O.Data)
                    S(k).Data = O.Data{k}.Data;
                    S(k).Idx = O.Data{k}.Idx;
                end
                O.Data = S;
            else
                error('%s: Unrecognized input data format. Must be a cell of struct with field Data and Idx.', mfilename);
            end
        end

        function reset(O)
            O.Ax = struct('Template', [], 'Confetti', [], 'Diversity', [], 'Entropy', []);
            O.ModData = O.Data;
            for j = 1:numel(O.Data)
                O.ModData(j).Group = O.Group(j);
            end
        end
        
        function saveData(O)
            if isempty(O.SaveDir)
                O.SaveDir = [];
            end
            if isempty(O.SaveDir)
                warning('%s: No save directory set', mfilename);
                return
            end
%             SaveHeader = ['Label', cellfun(@(x, y) sprintf('%d_Grp%d', x, y), num2cell(1:numel(O.Group)), num2cell(O.Group), 'un', 0)];
%             SaveData = [O.ModData.Label num2cell(O.ModData.Indiv)];
%             SaveName = fullfile(O.SaveDir, [O.DataName '.csv']);
%             writeDlmFile(vertcat(SaveHeader, SaveData), SaveName);
%             fprintf('%s: Save to file "%s".\n', mfilename, SaveName);
        end
        
        function savePlot(O)
            if isempty(O.SaveDir)
                warning('%s: No save directory set.', mfilename);
                return
            end
            if ~isempty(O.Gx) && O.isGxValid
                SaveName = fullfile(O.SaveDir, [O.DataName '.png']);
                savePlot(O.Gx, 'SaveAs', SaveName, 'DPI', 600);
            else
                warning('%s: No figure handle found.', mfilename);
                return
            end
        end
    end
    
    %Plotting methods
    methods
        
        function plotRarefaction(O)
        end
        
        function plotConfetti(O)
            Data = O.getDataAsCell;
            
            %Make sure the confetti map matrices are made
            if isempty(O.ConfettiMap) || numel(O.ConfettiMap) ~= numel(Data)
                O.ConfettiMap = cell(1, numel(Data));
            end
            for j = 1:numel(Data)
                if isempty(O.ConfettiMap{j}) || size(O.ConfettiMap{j}, 1) ~= size(O.ConfettiMap{j}, 2)
                    fprintf('%s: Making ConfettiMap for Data # %d.\n', mfilename, j)
                    O.ConfettiMap{j} = makeConfettiMap(Data{j});
                end
            end
            
            %plot the confetti maps in a MxN grid
            [~, ~, ~, GrpIdx] = unique2([O.ModData.Group]);
            R = max(cellfun('length', GrpIdx));
            C = numel(GrpIdx);
            
            if isempty(O.Ax.Confetti) || endsWith(class(O.Ax.Confetti(1)), 'GraphicsPlaceholder') || ~isvalid(O.Ax.Confetti(1))
                Gx = figure;
                O.Ax.Confetti = gobjects(R, C);
                j = 1;
                for r = 1:R
                    for c = 1:C
                        O.Ax.Confetti(r, c) = subplot(R, C, j);
                        j = j + 1;
                    end
                end
            else
                Gx = get(O.Ax.Confetti(1), 'parent');
            end
            
            %Place a black image on all as default
            for j = 1:numel(O.Ax.Confetti)
                plotConfetti(O.Ax.Confetti(j), 1, [0 0 0]);
            end
            %Draw the Confetti plots
            for j = 1:numel(GrpIdx)
                for k = 1:numel(GrpIdx{j})
                    plotConfetti(O.Ax.Confetti(k, j), O.ConfettiMap{GrpIdx{j}(k)});
                end
            end
            %Post modifcation
            resizeFigure(Gx, 2*C, 2*R);
            resizeSubplots(Gx, 'HorzSpacer', 0.01);
            centerFigureOnMonitor(Gx);
        end
        
        function saveConfetti(O, varargin)
            %Check if data was pooled by group
            IsPooled = numel(O.Data) ~= numel(O.ModData);
            PoolName = ternary(IsPooled, 'Pooled', 'File');
            
            %See if downsample was used
            NumelData = cellfun('length', {O.ModData.Data});
            IsDownSampled = all(NumelData == min(NumelData));
            DownSampleName = ternary(IsDownSampled, sprintf('.DownSample%d', min(NumelData)), '');
            
            FileName = sprintf('%s.Confetti.%s%s.png', O.DataName, PoolName, DownSampleName);
            SaveName = fullfile(O.SaveDir, FileName);

            try 
                Gx = get(O.Ax.Confetti(1), 'parent');
            catch
                O.plotConfetti;
                Gx = get(O.Ax.Confetti(1), 'parent');
            end
            savePlot(Gx, 'SaveAs', SaveName, varargin{:});
        end
        
        function saveIndivConfetti(O, varargin)
            %Check if data was pooled by group
            IsPooled = numel(O.Data) ~= numel(O.ModData);
            PoolName = ternary(IsPooled, 'Pooled', 'File');
            
            %See if downsample was used
            NumelData = cellfun('length', {O.ModData.Data});
            IsDownSampled = all(NumelData == min(NumelData));
            DownSampleName = ternary(IsDownSampled, sprintf('.DownSample%d', min(NumelData)), '');
            
            for j = 1:numel(O.ConfettiMap)
                if IsPooled
                    FileName = sprintf('%s.Confetti.Grp%d.%s%s.png', O.DataName, O.ModData(j).Group, PoolName, DownSampleName);
                else
                    FileName = sprintf('%s.Confetti.Grp%d.%s%d%s.png', O.DataName, O.ModData(j).Group, PoolName, j, DownSampleName);
                end
                SaveName = fullfile(O.SaveDir, FileName);
                Gx = plotConfetti(O.ConfettiMap{j});
                savePlot(Gx, 'SaveAs', SaveName, varargin{:});
                close(Gx);
            end
        end
        
        function plotCumDist(O)
            plotCDF(O)
        end           
        
        function saveCumDist(O, varargin)
            FileName = [O.DataName 'CDF.png'];
            SaveName = fullfile(O.SaveDir, FileName);
            try 
                Gx = get(O.Ax.CDF, 'parent');
            catch
                return
            end
            savePlot(Gx, 'SaveAs', SaveName, varargin{:});
        end
        
        function plotDiversity(O, varargin)
            if isempty(varargin)
                varargin = {'richness'};
            end
            Data = O.getDataAsCell;
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
            
            Mean = applyall('mean', Diversity, [O.ModData.Group]);
            Error = applyall('std', Diversity, [O.ModData.Group]);
            Ax = O.Ax.Diversity;
            if isempty(Ax) || ~isvalid(Ax)
                Ax = [];
            end
            Gx = plotErrorBox(Ax, Mean, 'Std', Error, 'IndivData', Diversity, 'Group', [O.ModData.Group]);
            Ax = get(Gx, 'children');
            O.Ax.Diversity = Ax;
            ylabel(Ax, 'Diversity');
            Ax.XTickLabel = varargin;
            resizeSubplots(Ax);
        end
        
        function saveTemplate(O, varargin)
            FileName = [O.DataName 'Scatter.png'];
            SaveName = fullfile(O.SaveDir, FileName);
            try 
                Gx = get(O.Ax.Template(1), 'parent');
            catch
                return
            end
            savePlot(Gx, 'SaveAs', SaveName, varargin{:});
        end
        
        function plotEntropy(O, varargin)
            if isempty(varargin)
                varargin = {'richness'};
            end
            Data = O.getDataAsCell;
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
            
            Mean = applyall('mean', Entropy, [O.ModData.Group]);
            Error = applyall('std', Entropy, [O.ModData.Group]);
            Ax = O.Ax.Simpson;
            if isempty(Ax) || ~isvalid(Ax)
                Ax = [];
            end
            Gx = plotErrorBox(Ax, Mean, 'Std', Error, 'IndivData', Entropy, 'Group', [O.ModData.Group]);
            Ax = get(Gx, 'children');
            O.Ax.Entropy = Ax;
            ylabel(Ax, 'Entropy');
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
        
        function plot(O, varargin)
            O.plotConfetti(varargin{:});
            O.plotTemplate(varargin{:});
            O.plotDiversity(varargin{:});
        end    
    end
  
    %set and get private properties
    methods       
        function set.Group(O, NewGroup)
            if nargin == 1 || isempty(NewGroup)
                O.PrivateGroup = repelem(1, 1, size(O.Data, 2)-1);
            elseif numel(NewGroup) == 1
                O.PrivateGroup = repelem(NewGroup, 1, size(O.Data, 2)-1);
            elseif numel(NewGroup) ~= size(O.Data, 2)
                error('%s: Number of new group elements differs from that of Data frequency columns.', mfilename);
            else
                O.PrivateGroup = NewGroup;
            end
            O.reset; %redo the grouping.
        end
        
        function Group = get.Group(O)
            Group = O.PrivateGroup;
        end 
    end
    
    methods (Access = private)
        function Data = getDataAsCell(O, Option)
            Data = cell(1, numel(O.ModData));
            if nargin < 2
                Option = '';
            end
            DoPercentage = endsWith(Option, 'perc', 'ignorecase', true);
            DoNormalize = startsWith(Option, {'norm', 'normalize'}, 'ignorecase', true);
            for f = 1:numel(O.ModData)
                if DoNormalize
                    Data{f} = O.ModData(f).Data ./ sum(O.ModData(f).Data);
                    if DoPercentage
                        Data{f} = Data{f} * 100;
                    end
                else
                    Data{f} = O.ModData(f).Data;
                end
            end
        end
    end
end