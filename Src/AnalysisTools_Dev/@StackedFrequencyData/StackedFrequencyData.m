%StackedFrequencyData is an object container for frequency plots for BRILIA. It
%handles the storage, manipulation, ploting, and saving of such data.
%
%  Methods
%    plot: plot the current working data set
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
%  NOTE: This is done PER FREQUENCY data of ALL REPERTOIRE
classdef StackedFrequencyData < DataInterface
    properties (Dependent)
        Group double    %grouping of data
    end
    
    properties (Access = public)
        XTickLabel cell %Cell array for X-axis
        CtrlGroup %Ctrl group
    end
     
    properties (GetAccess = public, SetAccess = private)
        GroupStat      struct %stores information about the groupwise test
        PrivateGroup   double %stores the current grouping of data, workaround dependency issues in MATLAB
        Width          double %width in inch of figure
        Height         double %height in inch of figure
    end
    
    methods (Access = public)
        function O = StackedData(varargin)
            P = inputParser;
            addRequired(P, 'Data', @(x) iscell(x) || isnumeric(x));
            addParameter(P, 'DataName', mfilename, @ischar);
            addParameter(P, 'Group', [], @isnumeric);
            addParameter(P, 'CtrlGroup', [], @isnumeric);
            addParameter(P, 'Width',  6, @(x) isnumeric(x) && x > 0.1);
            addParameter(P, 'Height', 3, @(x) isnumeric(x) && x > 0.1);
            parse(P, varargin{:});
            O.Width = P.Results.Width;
            O.Height = P.Results.Height;
            O.Data = P.Results.Data;
            O.DataName = P.Results.DataName;
            O.format; %Do this BEFORE assigning group, otherwise error can occur.
            O.Group = P.Results.Group;
            O.CtrlGroup = P.Results.CtrlGroup;
            O.reset; %Resets the ModData field
        end
        
        function TF = isValid(O)
            TF = iscell(O.Data) && ~isempty(O.Data);
        end
        
        function filter(O, Option)
            if strcmpi(Option, 'pool') %pool data by group
                disp('Pooling data by group.');
                [UnqGroup, ~, UnqIdx] = unique(O.Group);
                O.ModData.Indiv = splitapply(@sum, O.ModData.Indiv', UnqIdx)';
                O.ModData.Group = UnqGroup;
            end
        end
        
        function O = format(O)
            if any(size(O.Data)) == 1 && ~iscell(O.Data(1))%Data is Nx1 matrix or cell -> count immediately
                O.Data = countData(O.Data);
            elseif iscell(O.Data) 
                if size(O.Data, 2) >= 2 && all(all(cellfun('isclass', O.Data(:, 2:end), 'double'))) && numel(unique(cellfun('length', O.Data(:, 2:end)))) == 1 %Already conformed
                    O.Data = O.Data;
                elseif all(cellfun(@(x) min(size(x)), O.Data) == 1) %cell containing lists -> count per cell
                    CellData = cellfun(@countData, O.Data, 'un', 0);
                    O.Data = conformDist(CellData{:});
                elseif all(cellfun(@(x) size(x,2) == 2 && iscell(x), O.Data)) %cell contains lists of Mx2 cells
                    O.Data = conformDist(O.Data{:});
                
                end
            end
            if isempty(O.Data)
                error('%s: Data must be either:\n MxN cell array where Col1 = category, Col2+ = frequency\n 1xN cell array with each cell storing countable (Mx1) or counted (Mx2) data.', mfilename)
            end
        end              
        
        function O = normalize(O)
            O.ModData.Indiv = O.ModData.Indiv ./ sum(O.ModData.Indiv, 1);
        end
       
        function O = reset(O)
            O.ModData = struct('Indiv', [], 'Label', []);
            O.ModData.Indiv = cell2mat(O.Data(:, 2:end));
            O.ModData.Label = O.Data(:, 1);
            O.ModData.Group = O.Group;
        end
        
        function saveData(O)
            if isempty(O.SaveDir)
                O.SaveDir = [];
            end
            if isempty(O.SaveDir)
                warning('%s: No save directory set', mfilename);
                return
            end
            SaveHeader = {'BC', 'SC'}; 
            SaveData = num2cell(O.Data);
            SaveName = fullfile(O.SaveDir, [O.DataName '.csv']);
            writeDlmFile(vertcat(SaveHeader, SaveData), SaveName);
            fprintf('%s: Save to file "%s".\n', mfilename, SaveName);
        end
        
        function savePlot(O)
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
    end
    
    %Plotting methods
    methods
        function O = plot(O)
            if isempty(O.Ax) || ~(all(isvalid(O.Ax))) && numel(O.Ax) ~= 1
                figure;
                O.Ax = axes;
            end
            bar(O.Ax, O.ModData.Indiv', 'stacked');
            setAxes(O.Ax, 'XTickLabel',[O.ModData.Group], 'XTick', 1:numel([O.ModData.Group]))
            Label = strrep(O.ModData.Label, '_', ' ');
            Lx = flippedLegend(O.Ax, Label{:}, 'Location', 'eastoutside');
            O.color;
            O.addNumber;
            setPlotTickDecimal(O.Ax, -1, 1);
            resizeFigure(O.Gx, O.Width, O.Height);
            resizeSubplots(O.Gx);
            O.Ax.Position(3) = O.Ax.Position(3) - Lx.Position(3)*1.2;
        end        
        
        function O = plotPie(O)
            if isempty(O.Ax) || ~(all(isvalid(O.Ax)))
                figure;
                O.Ax = axes;
            end
            N = size(O.ModData.Indiv, 2);
            Px = cell(1, N);
            ColorMap = [0.01 0.01 0.01; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 1 1 1]; 
            for j = 1:N
                O.Ax(j) = subplot(1, N, j);
                Px{j} = pie(O.Ax(j), O.ModData.Indiv(:, j));
                colormap(O.Ax(j), ColorMap);
                ColorIdx = round(linspace(1, size(ColorMap, 1), size(O.ModData.Indiv, 1)));
                for k = 1:2:numel(Px{j})
                    %Object slice
                    Ox = Px{j}(k);
                    Ox.ZData = -ones(size(Ox.XData)); %So that it's behind text
                    Ox.LineWidth = 1.5;
                    %Text
                    Tx = Px{j}(k+1);
                    Tx.FontName = 'Arial';
                    Tx.FontSize = 12;
                    Tx.FontWeight = 'Bold';
                    Tx.Position = Tx.Position * 0.5;
                    Color = 1 - ColorMap(ColorIdx((k+1)/2), :);
                    Color(Color > 1) = 1;
                    Tx.Color = Color;
                    Tx.String = strrep(Tx.String, '%', '');
                    title(O.Ax(j), sprintf('%d', O.ModData.Group(j)), 'FontName', 'Arial', 'FontSize', 24, 'FontWeight', 'bold');
                end
            end
            resizeFigure(gcf, O.Width, O.Height);
            resizeSubplots(gcf);
        end
        
        
        function savePie(O)
            SaveName = fullfile(O.SaveDir, [O.DataName '.pie.png']);
            Gx = get(O.Ax(1), 'parent');
            savePlot(Gx, 'SaveAs', SaveName);
        end
            
        function color(O, ClrMap)
            if nargin == 1
                ClrMap = getStdColor(1.4);
                ClrMap(1, :) = []; %Black is unappealing
            end
            
            if ~isempty(O.Ax) && isvalid(O.Ax)
                Bx = findobj(O.Ax, 'Type', 'Bar');
            else
                return
            end
            
            c = 1;
            for j = 1:numel(Bx)
                Bx(j).FaceColor = ClrMap(j, :);
                c = incr(c, size(ClrMap, 1), 1);
            end
        end
        
        function setLegend(O, varargin)
            if ~O.isGxValid; return; end
            Legend = findobj(O.Gx, 'Type', 'Legend');
            set(Legend, varargin{:});
        end
        
        function addNumber(O, varargin)
            if ~O.isGxValid; return; end
            O.delNumber;
            if nargin > 1
                Fmt = varargin{1};
            else
                Fmt = '%0.2f';
            end
            Y = O.ModData.Indiv;
            X = repmat(1:size(O.ModData.Indiv, 2), size(Y, 1), 1);
            Ypos = cumsum(Y, 1) - Y/2;
            for k = 1:numel(Y)
                text(O.Ax, X(k), Ypos(k), sprintf(Fmt, Y(k)), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle', 'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'bold');
            end
        end
        
        function delNumber(O)
            if ~O.isGxValid; return; end
            Tx = findobj(O.Gx, 'type', 'text');
            delete(Tx)
        end
    end
  
    %set and get private properties
    methods       
        function O = set.Group(O, NewGroup)
            if nargin == 1 || isempty(NewGroup)
                O.PrivateGroup = repelem(1, 1, size(O.Data, 2)-1);
            elseif numel(NewGroup) == 1
                O.PrivateGroup = repelem(NewGroup, 1, size(O.Data, 2)-1);
            elseif numel(NewGroup) ~= size(O.Data, 2)-1
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
end