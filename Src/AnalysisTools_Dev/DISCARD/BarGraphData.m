%BarGraphData is an object container for frequency plots for BRILIA. It
%handles the storage, manipulation, ploting, and saving of such data.
%
%  Data Input
%      - Px2 array where N(:,1) are values and N(:,2) are frequencies
%      - Qx1 array where each element is a unique observation
%
%CODING_NOTE: 
classdef BarGraphData < DataInterface
    properties (Access = public)
        XTickLabel cell %Cell array for X-axis
    end
     
    properties (GetAccess = public, SetAccess = private)
        MeanHandle     matlab.graphics.Graphics %Mean handles
        ErrorHandle    matlab.graphics.Graphics %Mrror handles
        SigTestHandle  matlab.graphics.Graphics %Text handles of stat. sig. pairwise comparisons
        GroupStat      struct %Results of groupwiseTest
        Alpha          double %Alpha value for doing groupwise comparisons
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
    
    methods (Access = public)
        function O = BarGraphData(varargin)
            %Creates a BarGraphData object
            O = O@DataInterface; %Use superclass constructor (REQUIRED)
            O.setGx('BarGraph', []);
            
            %Extract and add the Data to object
            DataIdx = find(cellfun(@ischar, varargin), 1)-1;
            if isempty(DataIdx); DataIdx = numel(varargin); end
            O.Data = BarGraphData.format(varargin{1:DataIdx});
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
            %Adds a data point to BarGraphData
            Data = BarGraphData.format(Data);
            O.add@DataInterface(Data, varargin{:});
            O.reset;
        end
        
        function del(O, Number)
            %Deletes a data point from BarGraphData
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
            %Resets the BarGraphData's internal data structure
            
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
            O.ModData.Group = applyall('mean', O.ModData.Indiv, O.Group);
            O.ModData.Error = applyall('std',  O.ModData.Indiv, O.Group);
            O.ModData.PairH = [];
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
    
    %Subclass methods
    methods
        function rescale(O)
            %Rescales the frequencies to match the mean of the CtrlGroup values
            Idx = find(O.Group == O.CtrlGroup);
            Idx = ternary(isempty(Idx), 1, Idx);
            [O.ModData.Indiv, O.ModData.MultFactor] = rescaleDist(O.ModData.Indiv, O.ModData.Indiv(:, Idx));
            O.ModData.Group = applyall('mean', O.ModData.Indiv, O.Group);
            O.ModData.Error = applyall('std',  O.ModData.Indiv, O.Group);
            O.ModData.PairH = [];
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
            O.ModData.Group = applyall('mean', O.ModData.Indiv, O.Group);
            O.ModData.Error = applyall('std',  O.ModData.Indiv, O.Group);
            O.ModData.PairH = [];
        end

        function normalize(O)
            %Normalizes frequencies so that sum of frequencies = 1
            O.ModData.Indiv = O.ModData.Indiv./sum(O.ModData.Indiv);
            O.ModData.Group = applyall('mean', O.ModData.Indiv, O.Group);
            O.ModData.Error = applyall('std',  O.ModData.Indiv, O.Group);
            O.ModData.PairH = [];
            O.ModData.MultFactor = ones(1, size(O.Data, 2));
        end
        
        function setError(O, Option)
            %Sets the error bars according to "minmax" of data or "std" around mean of data.
            switch lower(Option)
                case 'minmax'
                    MaxData = applyall('max', O.ModData.Indiv, O.Group);
                    MinData = applyall('min', O.ModData.Indiv, O.Group);
                    MidData = (MaxData + MinData) / 2;
                    O.ModData.Error = MaxData - MidData;
                    O.ModData.Group = MidData;
                case 'std'
                    O.ModData.Group = applyall('mean', O.ModData.Indiv, O.Group);
                    O.ModData.Error = applyall('std',  O.ModData.Indiv, O.Group);
            end
        end
        
        %Get the significant difference between group comparisons
        function groupwiseTest(O, varargin)
            %Conducts a group-wise comparison of frequencies
            [~, Gall]= groupwiseTest(O.ModData.Indiv, O.Group, 'Alpha', O.Alpha, varargin{:});
            O.ModData.PairH = Gall(:, 3:end)' <= O.Alpha;
        end
    end
    
    %Plot methods
    methods
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
        
        function plot(O, TopN)
            %Plots all data or the top N bars
            if ~isfield(O.ModData, 'Indiv')
                disp('There is no data left to plot.');
                return
            end
            if nargin > 1 && size(O.ModData.Indiv, 1) >= TopN
                DataSum = max(O.ModData.Indiv, [], 2);
                [~, Idx] = sort(DataSum, 'descend');
                ShowIdx = sort(Idx(1:TopN));
            else
                ShowIdx = [];
            end
            if ~isfield(O.ModData, 'PairH')
                O.ModData.PairH = [];
            end
            if ~isnumeric(O.Gx.BarGraph) && ~isvalid(O.Gx.BarGraph)
                O.Gx.BarGraph = [];
            end
            [O.Gx.BarGraph, O.MeanHandle, O.ErrorHandle, O.SigTestHandle] = plotErrorBox(O.Gx.BarGraph, O.ModData.Group, 'STD', O.ModData.Error, 'PairH', O.ModData.PairH, 'XTickLabel', O.ModData.Label, 'IndivData', O.ModData.Indiv, 'Group', O.Group, 'ShowIdx', ShowIdx);
            if ~isempty(O.DataName)
                TitleName = strrep(sprintf('%s', strrep(O.DataName, '_', ' ')), '_', '\_');
                setAxes(O.Gx.BarGraph, 'XLabel-String', O.DataName, 'Title-String', TitleName);
            end
            resizeFigure(O.Gx.BarGraph, O.Width, O.Height);
            resizeSubplots(O.Gx.BarGraph);
        end

        function colorBars(O, RGB)
            %Colors each main bar a certain color
            arrayfun(@(x) set(x, 'FaceColor', RGB), O.MeanHandle);
        end
 
        function colorErrorBars(O, RGB)
            %Colors each error bar a certain color
            arrayfun(@(x) set(x, 'FaceColor', RGB), O.ErrorHandle);
        end
        
        function showErrorBars(O)
            %Shows error bars
            arrayfun(@(x) set(x, 'Visible', 'on'), O.ErrorHandle);
        end

        function hideErrorBars(O)
            %Hides error bars
            arrayfun(@(x) set(x, 'Visible', 'off'), O.ErrorHandle);
        end

        function showBars(O)
            %Shows main bars
            arrayfun(@(x) set(x, 'Visible', 'on'), O.MeanHandle);
        end
        
        function hideBars(O)
            %Hides main bars
            arrayfun(@(x) set(x, 'Visible', 'off'), O.MeanHandle);
        end
    end
end