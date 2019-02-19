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
    methods (Access = public)
        function O = ConvergenceData(varargin)
            O = O@DataInterface; %Use superclass constructor (REQUIRED)
            
            %Extract and add the Data to object
            DataIdx = find(cellfun(@ischar, varargin), 1)-1;
            if isempty(DataIdx); DataIdx = numel(varargin); end
            O.Data = ConvergenceData.format(varargin{1:DataIdx});
            varargin = varargin(DataIdx+1:end);
            
            %Parse the remaining inputs
            P = inputParser;
            P.addParameter('DataName', mfilename, @ischar);
            P.addParameter('Group', [], @isnumeric);
            P.addParameter('GroupName', {}, @iscell);
            P.addParameter('CtrlGroup',  1, @isnumeric);
            P.addParameter('Width',  6, @(x) isnumeric(x) && x > 1);
            P.addParameter('Height', 6, @(x) isnumeric(x) && x > 1);
            P.parse(varargin{:});
            
            %Update the object based on parsed inputs
            O.DataName = P.Results.DataName;
            O.Group = P.Results.Group;
            O.GroupName = P.Results.GroupName;
            O.CtrlGroup = P.Results.CtrlGroup;
            O.Width  = P.Results.Width;
            O.Height = P.Results.Height;
            
            O.reset;  %Resets the ModData field
        end
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
                %Output data should be a Mx2 cell array, where col1 are unique char ID, and col2 are weights
                %NOTE: Col1 can have duplicates, for example, to indicate different clonotypes with same CDR3
                if isstruct(Data) && isfield(Data, 'Data') %Unwrap Data from stucture only
                    Data = Data.Data;
                    return
                end
                if ~(iscell(Data) && size(Data, 2) == 2 && all(cellfun(@(x, y) ischar(x) && isnumeric(y), Data(:, 1), Data(:, 2))))
                    error('%s: Unrecognized input data format. Must be a Mx2 cell or a struct with Data field. Col1 must be char, and Col2 must be numbers.', mfilename);
                end
            end
        end
    end
    
    methods
        
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
        
        function pool(O)
        end

        function open(O)
        end
                
        function reset(O)
            O.ModData = struct;
            O.ModData.Data = {O.Data.Data};
        end

        function save(O, FMT)
            if isempty(O.SaveDir)
                O.SaveDir = [];
            end
            if isempty(O.SaveDir)
                warning('%s: No save directory set', mfilename);
                return
            end
            
            if nargin < 2
                FMT = '%0.3f';
            end
            
            if ~isfield(O.ModData, 'Heatmap'); return; end
            GrpAndFile = cellfun(@(x, y) sprintf('%s (%d)', x, y), O.GroupName(:)', num2cell(1:numel(O.Data)), 'un', 0);
            SaveData = cell(numel(O.Data) + 1);
            SaveData(1, 2:end) = GrpAndFile;
            SaveData(2:end, 1) = GrpAndFile;
            
            Fields = fieldnames(O.ModData.Heatmap);
            for f = 1:numel(Fields)
                Unit = ternary(contains(Fields{f}, 'Overlap', 'ignorecase', 1), '%', 'Unitless');
                Multiplier = ternary(strcmpi(Unit, '%'), 100, 1);
                SaveData{1} = sprintf('Group (Subj#) [%s]', Unit);
                SaveName = fullfile(O.SaveDir, sprintf('%s.%s.csv', O.DataName, Fields{f}));
                NumData = num2cell(O.ModData.Heatmap.(Fields{f}) * Multiplier);
                StrNumData = cellfun(@(x) sprintf(FMT, x), NumData, 'un', 0);
                SaveData(2:end, 2:end) = StrNumData;
                writeDlmFile(SaveData, SaveName);
                fprintf('Saved file: "%s".\n', mfilename, SaveName);
            end
        end
        
        function savePlot(O)
            %Saves all plot to the path "SaveDir\DataName.png"
            if isempty(O.SaveDir)
                warning('%s: No save directory set.', mfilename);
                return
            end
            Fields = fieldnames(O.Gx);
            for f = 1:numel(Fields)
                if isempty(O.Gx.(Fields{f})) || ~isvalid(O.Gx.(Fields{f}))
                    continue
                else
                    SaveName = fullfile(O.SaveDir, sprintf('%s.%s.png', O.DataName, Fields{f}));
                    savePlot(O.Gx.(Fields{f}), 'SaveAs', SaveName, 'DPI', 600);
                    fprintf('Saved file: "%s"\n', SaveName);
                end
            end
        end
    end
    
    %Plotting methods
    methods
        function plot(O, Fields)
            if ~isfield(O.ModData, 'Heatmap') || isempty(O.ModData.Heatmap) || numel(fieldnames(O.ModData.Heatmap)) == 0
               O.getHeatmapData;
            end
            if nargin < 2
                Fields = fieldnames(O.ModData.Heatmap);
            end
            
            for f = 1:numel(Fields)
                FigName = Fields{f};
                if ~O.isGxValid(FigName)
                    O.setGx(FigName, figure);
                end
                FMT = ternary(contains(Fields{f}, 'overlap', 'ignorecase', true), '%0.1f', '%0.2f');
                Data = O.ModData.Heatmap.(Fields{f});
                Data(1:size(Data,1)+1:end) = NaN;
                heatmap2(Data, O.Group, O.Group, FMT, 'Colormap', cool, 'Parent', O.Ax.(FigName), 'ColorBar', [], 'NaNColor', [0 0 0], 'ShowAllTicks', true);
            end
        end
    end
     
    methods
        %HeatmapData is a MxM matrix of overlap 0.2f % between repertoire
        %CDR3 sequences. The fraction is always in M overlatp betwen Rth
        %row and Cth column repertoire.
        %EX:  Heatmap = [100  20;  <- 20% of Rep 2 overlaps with Rep 1
        %                 30 100]  <- 30% of Rep 1 overlaps with Rep 2
        function getHeatmapData(O)
            O.ModData.Heatmap = calcConvMatrix(O.ModData.Data{:});
        end
    end
end