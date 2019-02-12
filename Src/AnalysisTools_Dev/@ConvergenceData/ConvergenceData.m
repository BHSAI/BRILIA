%ConvergenceData is an object container for many different plots of
%identifying similary between repertoires. It handles the storage,
%manipulation, ploting, and saving of such data.
%
%  Methods
%    plot: plot ALL types of convergence plots for the current working data set
%    plotCDR3
%    plotKmer(k) where k is the k-mer
%    plotCDR3Prop(AtoP) for reduced-letter, where AtoP is a 21x1 letter code for conversion. default is by convAA2Prop
%    plotKmerProp(AtoP) for reduced-letter kmers
%    plotDendrogram plots reperotire similarity by ...
%    plotCDR3LineConvergence 
%    
%    savePlot: saves the current plot
%    saveData: saves the current data
%
%  Properties
%    Ax: struct of axes
%    Gx: struct of figures
%    ModData.CDR3: MxM square matrix of overlap CDR3
%    ModData.Kmer: MxMxY matrix of overlap Kmer Freq for Y unique 
%    ModData.CDR3Prop: MxM 
%    ModData.KmerProp: MxMxY matrix of overlap kmer properties for Y unique ...
%
%  NOTE: this is done for ALL lineages PER REPERTOIRE
%
classdef ConvergenceData < DataInterface
    properties (Dependent)
        Group double    %grouping of data
    end
    
    properties (Access = public)
        CtrlGroup %Ctrl group
        Level %Clustering level of template: Clone, Clonotype, CDR3, Kmer, etc
        Color %Color map Mx3 matrix
        Width          double %width in inch of figure
        Height         double %height in inch of figure
        Heatmap
    end
     
    properties (GetAccess = public, SetAccess = private)
        GroupStat      struct %stores information about the groupwise test
        PrivateGroup   double %stores the current grouping of data, workaround dependency issues in MATLAB
    end
    
    methods (Access = public)
        function O = ConvergenceData(varargin)
            P = inputParser;
            addRequired(P, 'Data', @(x) iscell(x) || isnumeric(x));
            addParameter(P, 'DataName', mfilename, @ischar);
            addParameter(P, 'Group', [], @isnumeric);
            addParameter(P, 'CtrlGroup', [], @isnumeric);
            addParameter(P, 'Width',  6, @(x) isnumeric(x) && x > 0.1);
            addParameter(P, 'Height', 6, @(x) isnumeric(x) && x > 0.1);
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
        
        function TF = isValid(O)
            TF = iscell(O.Data) && ~isempty(O.Data);
        end
        
        function filter(O, varargin)
            P = inputParser;
            P.addParameter('Level', 'Clone', @(x) ismember(lower(x), {'clone', 'ac', 'bc', 'sc'}));
            P.addParameter('Alphabet', 'AA', @(x) ismember(lower(x), {'aa', 'prop'}));
            P.addParameter('GroupBy', '', @ischar);
            parse(P, varargin{:});
            O.Level = P.Results.Level;
            Alphabet = P.Results.Alphabet;
            GroupBy = P.Results.GroupBy;
            O.reset; 
            
            for f = 1:numel(O.Data)
                G = getGrpIdx(O.Data(f).GrpNum, O.Data(f).Template, O.Level);
                O.ModData(f).GrpNum = [G.GrpNum]';
                O.ModData(f).Template = [G.Template]';
                Fields = fieldnames(O.Data);
                Fields = Fields(endsWith(Fields, 'CDR3', 'ignorecase', true));
                if strcmpi(Alphabet, 'prop')
                    for j = 1:numel(Fields)
                        O.ModData(f).(Fields{j}) = cellfun(@convAA2PropMEX, O.ModData(f).(Fields{j})([G.Idx]'), 'un', 0);
                    end
                else
                    for j = 1:numel(Fields)
                        O.ModData(f).(Fields{j}) = O.Data(f).(Fields{j})([G.Idx]');
                    end
                end
                if isfield(O.ModData, GroupBy) %Add up the templates according to the field used : GrpNum or xCDR3
                    if strcmpi(GroupBy, 'GrpNum')
                        [O.ModData(f).Template, Idx] = regroupTemplate(O.ModData(f).Template, O.ModData(f).GrpNum);
                        Idx1 = cellfun(@(x) x(1), Idx);
                        O.ModData(f).GrpNum = O.ModData(f).GrpNum(Idx1);
                        for j = 1:numel(Fields)
                            O.ModData(f).(Fields{f}) = O.ModData(f).(Fields{j})(Idx1);
                        end                        
                    elseif strcmpi(GroupBy, 'hCDR3')
                        [O.ModData(f).Template, Idx] = regroupTemplate(O.ModData(f).Template, O.ModData(f).hCDR3);
                        Idx1 = cellfun(@(x) x(1), Idx);
                        O.ModData(f).GrpNum = O.ModData(f).GrpNum(Idx1);
                        for j = 1:numel(Fields)
                            O.ModData(f).(Fields{f}) = O.ModData(f).(Fields{j})(Idx1);
                        end                        
                    end
                end
            end
        end
        
        function format(O)
            %All data should be of Mx1 cell of struct of .Data, .Idx. This will format everything into a single nonscalar struct.
            if isstruct(O.Data) && isfield(O.Data, 'Data') && isfield(O.Data, 'Idx')
                return
            elseif iscell(O.Data) && all(cellfun('isclass', O.Data, 'struct'))
                Fields = fieldnames(O.Data{1});
                S(1:numel(O.Data)) = struct;
                for k = 1:numel(O.Data)
                    for f = 1:numel(Fields)
                        S(k).(Fields{f}) = O.Data{k}.(Fields{f});
                    end
                end
                O.Data = S;
            else
                error('%s: Unrecognized input data format. Must be a cell of struct with field Data and Idx.', mfilename);
            end
        end

        function reset(O)
            O.ModData = O.Data;
            O.Heatmap = struct;
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
       
        function plotHeatmap(O, Field)
            if isempty(O.Heatmap) || numel(fieldnames(O.Heatmap)) == 0
               O.getHeatmapData;
            end
            if nargin < 2
                Field = fieldnames(O.Heatmap);
            end
            
            if ~isfield(O.Gx, 'Heatmap') || isempty(O.Gx.Heatmap) || ~isvalid(O.Gx.Heatmap)
                O.Gx.Heatmap = figure;
            end
            figure(O.Gx.Heatmap)
            
            %Initialize Axes
            for f = 1:numel(Field)
                if contains(Field{f}, 'ovlp', 'ignorecase', true)
                    Data = O.Heatmap.(Field{f}) * 100; %Want % for overlap %
                    Fmt = '%0.1f';
                    AppendTitle = ' (%)';
                    ColorMap = cool;
                else
                    Data = O.Heatmap.(Field{f});
                    Fmt = '%0.2f';
                    AppendTitle = '';
                    ColorMap = cool;
                    ColorMap = 1.2 - ColorMap(:, [2 3 1]);
                    ColorMap(:, 2) = ColorMap(:, 2)*1.3;
                    ColorMap(:, 3) = ColorMap(:, 3)*2;
                    ColorMap(ColorMap>1) = 1;
                end
                Data(1:size(Data,1)+1:end) = NaN;
                
                
%                 Ax = subplot(2, 2, f);
%                 heatmap2(Data, O.Group, O.Group, Fmt, 'Parent', Ax, 'ColorMap', ColorMap, 'ColorBar', [], 'NaNColor', [0 0 0], 'ShowAllTicks', true);
%                 title(Ax, [O.DataName '-' Field{f} AppendTitle]);
                
                %Temporary save indivi figures.
                FileName = [O.DataName '.' Field{f} '.Heatmap.png'];
                SaveName = fullfile(O.SaveDir, FileName);
                TGX = figure;
                TAX = axes;
                heatmap2(Data, O.Group, O.Group, Fmt, 'Parent', TAX, 'ColorMap', ColorMap, 'ColorBar', [], 'NaNColor', [0 0 0], 'ShowAllTicks', true);
                title(TAX, [O.DataName ' - ' Field{f} AppendTitle]);
                resizeFigure(TGX, 5, 5);
                setAxes(TAX, 'FontName', 'Arial', 'FontSize', 16);
                resizeSubplots(TAX);
                centerFigureOnMonitor(TGX)

                savePlot(TGX, 'SaveAs', SaveName)
            end
            
%             pause(0.5)
%             O.resizeHeatmapFigure(O.Width, O.Height)
%             pause(0.5)
%             O.resizeHeatmapSubplot('FigSpacer', 0.01, 'HorzSpacer', 0.01', 'VertSpacer', 0.01', 'ScaleVertical', 'y');
%             centerFigureOnMonitor(O.Gx.Heatmap); %Needed to get the text to show for somereason
%             
%             drawnow
            %CODING_NOTE: heatmap does not yet have the Hx.XDisplayText option. Hence, the follow heatmap2 usage.
        end
        
        function resizeHeatmapSubplot(O, varargin)
            if isfield(O.Gx, 'Heatmap') && ~isempty(O.Gx.Heatmap) && isvalid(O.Gx.Heatmap)
                varargin = cleanCommandLineInput(varargin{:});
                resizeSubplots(O.Gx.Heatmap, varargin{:});
            end
        end
        
        function resizeHeatmapFigure(O, varargin)
            if isfield(O.Gx, 'Heatmap') && ~isempty(O.Gx.Heatmap) && isvalid(O.Gx.Heatmap)
                varargin = cleanCommandLineInput(varargin{:});
                resizeFigure(O.Gx.Heatmap, varargin{:});
            end
        end
        
        function saveHeatmap(O, varargin)
            FileName = [O.DataName '.Heatmap.png'];
            SaveName = fullfile(O.SaveDir, FileName);
            try 
                Gx = O.Gx.Heatmap;
            catch
                return
            end
            savePlot(Gx, 'SaveAs', SaveName, varargin{:});
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
    
    methods (Access = public)    
        %HeatmapData is a MxM matrix of overlap 0.2f % between repertoire
        %CDR3 sequences. The fraction is always in M overlatp betwen Rth
        %row and Cth column repertoire.
        %EX:  Heatmap = [100  20;  <- 20% of Rep 2 overlaps with Rep 1
        %                 30 100]  <- 30% of Rep 1 overlaps with Rep 2
        function getHeatmapData(O)
            InData = arrayfun(@(x) [x.hCDR3 num2cell(x.Template)], O.ModData, 'un', 0);
            O.Heatmap = calcConvMatrix(InData{:});
        end

    end
end