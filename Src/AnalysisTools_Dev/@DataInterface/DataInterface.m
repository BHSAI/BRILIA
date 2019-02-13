%DataInterface ensures that all Data types enforce contains basic, common
%functions for handling, plotting, and storing data.

classdef DataInterface < handle   
    
    %Main parameters of the Data objects used for plotting and analysis
    properties (Access = public, SetObservable, AbortSet)
        Data             %Raw data stored as a non-scalar structure. Fields are Data, Group, GroupName.
        ModData   struct %Modified data that is used to make the plots. Subclasses will have their own implementation.
        DataName  char   %File name to use for saving files
        SaveDir   char   %Folder name to use for saving files
        Default   struct %Default structure for plotting data, such as FontName, FontSize
        Width     double %Default figure width
        Height    double %Default figure height
    end
       
    %Dependent variables for allowing multiple private variables to be
    %modified (MATLAB workaround)
    properties (Access = public, Dependent)
        Ax               %Axes of handle
        Gx               %Figure handle
        Group double     %Group number for each data
        GroupName cell   %Group name of data of same size as Group
        CtrlGroup double %Control gorup number
    end
    
    %Private variables for facilitating value storage in MATLAB objects
    properties (Access = protected, SetObservable, AbortSet)
        PrivAx struct = struct('Ax', []) %Private axes structure of handles
        PrivGx struct = struct('Gx', []) %Private figure structure of handles
        PrivCtrlGroup double %Group number of the "control" group
    end
     
    %Common methods that each subclass must implement
    methods (Abstract)
        %Filters data to plot and updates ModData
        filter(O, Option)
        %Formats the data correctly, or returns error is not possible (STATIC)
        format(O, varargin)
        %Plots the data
        plot(O, varargin)
        %Pools data according to the assigned Group
        pool(O)
        %Resets the data back to raw data
        reset(O)
        %Opens a JSON file containing the data
        open(O, varargin)
        %Saves the data in JSON format
        save(O, varargin)
        %Saves the plot
        savePlot(O, varargin)
    end
    
    %Superclass methods that should be used before subclass overrides it
    methods 
        function O = DataInterface
        %Constructs a Data object with default setting and listeners.
        %
        %CODING_NOTE: All subclass should summon this constructor via
        %O = O@DataInterface, where "O" is the subclass object
            O.Data = struct('Data', [], 'Group', 1, 'GroupName', '1');
            O.ModData = struct;
            O.DataName = 'DataObject';
            O.SaveDir = '';
            O.Default = struct('FontName', 'Arial', 'FontSize', 12);
            O.Width = 5;
            O.Height = 5;
            O.CtrlGroup = 1;
            O.addlistener('PrivAx', 'PostSet', @O.handleChangedAx);
            O.addlistener('PrivGx', 'PostSet', @O.handleChangedGx);
        end
        
        function add(O, Data, varargin)
            %Adds a data point to a data set 
            if nargin < 3
                GroupOrName = 1;
            else
                GroupOrName = varargin{1};
            end
            if isnumeric(GroupOrName) %Group number was provided. Find Group name.
                Group = GroupOrName;
                GroupIdx = find(Group == O.Group, 1);
                GroupName = ternary(isempty(GroupIdx), num2str(GroupOrName), O.GroupName(GroupIdx)); %#ok<*PROPLC>
            else %Group name was provided. Find Group number.
                GroupName = GroupOrName;
                GroupIdx = find(strcmpi(GroupName, O.GroupName), 1);
                Group = ternary(isempty(GroupIdx), max(O.Group)+1, O.Group(GroupIdx));
            end
            k = numel(O.Data) + 1;
            O.Data(k).Data = Data{1};
            O.Data(k).Group = Group; 
            O.Data(k).GroupName = GroupName;        
        end
        
        function del(O, Number)
            %Deletes a data point defined by the Nth index
            if islogical(Number)
                Number = find(Number);
            end
            if max(Number) > numel(O.Data)
                warning('%s: The deletion number cannot exceed the number of data. No action performed.', mfilename);
                return
            elseif min(Number) < 1
                warning('%s: The deletion number cannot be lower than 1. No action performed.', mfilename);                
                return
            else
                Number = unique(Number);
            end
            if numel(Number) == numel(O.Data) %Special case, reset O.Data. Cannot delete everything.
                O.Data = struct('Data', [], 'Group', 1, 'GroupName', 1);
            else
                O.Data(Number) = [];
            end
        end
    end
    
    %Set and Get codes
    methods
        
        function IsValid = isAxValid(O, Field) %Is Gx valid?
            IsValid = isfield(O.Ax, Field) && ~isempty(O.Ax.(Field)) && isvalid(O.Ax.(Field));
        end
        function setAx(O, Field, NewAx)
            %Sets the axes handle of a particular field            
            if isa(NewAx, 'matlab.graphics.axis.Axes') || isempty(NewAx)
                O.PrivAx.(Field) = NewAx;
            else
                warning('%s: Input must be a MATLAB axes handle but is a %s instead. No changes made.', mfilename, class(NewAx));
            end
            O.handleChangedAx;
        end       
        function set.Ax(O, NewAx)
            if isstruct(NewAx) && all(structfun(@(x) isa(x, 'matlab.graphics.axis.Axes') || isempty(x), NewAx))
                O.PrivAx = NewAx;
            elseif isa(NewAx, 'matlab.graphics.axis.Axes')
                O.PrivAx.Ax = NewAx;
            else
                error('%s: Invalid input. Expected an axes handle or structure of axes handles.', mfilename);
            end
        end
        function Ax = get.Ax(O)
            Ax = O.PrivAx;
        end
        function Ax = getAx(O, Field) 
            %Gets the axes handle of a particular field 
            if isfield(O.PrivAx, Field)
                Ax = O.PrivAx.(Field);
            else
                Ax = [];
                warning('%s: Could not find field of Ax, "%s".', mfilename, Field);
            end
        end        
        
        function IsValid = isGxValid(O, Field) %Is Gx valid?
            IsValid = isfield(O.Gx, Field) && ~isempty(O.Gx.(Field)) && isvalid(O.Gx.(Field));
        end
        function setGx(O, Field, NewGx)
            %Sets the figure handle of a particular field          
            if isa(NewGx, 'matlab.ui.Figure') || isempty(NewGx)
                O.PrivGx.(Field) = NewGx;
            else
                warning('%s: Input must be a MATLAB figure handle but is a %s instead. No changes made.', mfilename, class(NewGx));
            end
            O.handleChangedGx;
        end               
        function set.Gx(O, NewGx)
            if isstruct(NewGx) && all(structfun(@(x) isa(x, 'matlab.ui.Figure') || isempty(x), NewGx))
                O.PrivGx = NewGx;
            elseif isa(NewGx, 'matlab.ui.Figure')
                O.PrivGx.Gx = NewGx;
            else
                error('%s: Invalid input. Expected a figure handle or structure of figure handles.', mfilename);
            end
        end
        function Gx = get.Gx(O)
            Gx = O.PrivGx;
        end       
        function Gx = getGx(O, Field)
            %Gets the figure handle of a particular field            
            if isfield(O.PrivGx, Field)
                Gx = O.PrivGx.(Field);
            else
                Gx = [];
                warning('%s: Could not find field of Gx, "%s".', mfilename, Field);
            end
        end
       
        function set.CtrlGroup(O, CtrlGroup)
            if isempty(CtrlGroup)
                CtrlGroup = 1;            
            else
                if isinf(CtrlGroup) 
                    CtrlGroup = numel(O.Data);
                elseif CtrlGroup < 1
                    CtrlGroup = 1;
                else
                    CtrlGroup = round(CtrlGroup(1));
                end
            end
            O.PrivCtrlGroup = CtrlGroup;
        end
        function CtrlGroup = get.CtrlGroup(O)
            CtrlGroup = O.PrivCtrlGroup;
        end
        
        function set.SaveDir(O, NewDir)
            if nargin < 2 || isempty(NewDir)
                NewDir = uiputdir2;
            end
            if iscell(NewDir) %Only get the first one
                O.SaveDir = NewDir{1};
            else
                O.SaveDir = NewDir;
            end
            %CODING_NOTE: Let subclasses make the output dir when needed.
            %Do NOT use mkdir here unless some files will actually go in.
        end        
        
        function set.Data(O, Data)
            if isempty(Data)
                return
            elseif isstruct(Data) %Probabily initializing Data
                O.Data = Data;
                return 
            end
            IsCellOfCell = iscell(Data) && all(cellfun(@iscell, Data));
            if IsCellOfCell
                for j = 1:numel(Data)
                    O.Data(j).Data = Data{j};
                end
            else
                O.Data.Data = Data;
            end
        end
        
        function set.DataName(O, Name)
            if ~ischar(Name)
                warning('%s: DataName input must be a char. No changes made.', mfilename);
            else
                O.DataName = Name;
            end
        end
        
        function set.Group(O, Group)
            if numel(Group) == 0
                Group = repelem(1, 1, numel(O.Data));
            elseif numel(Group) == 1
                Group = repelem(Group, 1, numel(O.Data));
            elseif numel(Group) ~= numel(O.Data)
                error('%s: Number of Group differs from that of Data.', mfilename);
            end
            GroupC = num2cell(Group);
            [O.Data.Group] = deal(GroupC{:});
        end
        function Group = get.Group(O)
            Group = [O.Data.Group]';
        end   
        
        function set.GroupName(O, GroupName)
            if isempty(GroupName)
                GroupName = [O.Data.Group];
            end
            if isnumeric(GroupName)
                GroupName = cellfun(@num2str, num2cell(GroupName), 'un', false);
            end
            if numel(GroupName) ~= numel(O.Data)
                error('%s: Cannot have unequal number GroupName and Data.', mfilename)
            end
            [O.Data.GroupName] = deal(GroupName{:});
            %See if the Group and GroupName uniqueness matches
            [UnqGroupName, ~, UnqGroupIdx] = unique(GroupName);
            UnqGroup = unique(O.Group);
            if numel(UnqGroupName) ~= numel(UnqGroup) %Update Group
                UnqGroupIdxC = num2cell(UnqGroupIdx);
                [O.Data.Group] = deal(UnqGroupIdxC{:});
            end
        end       
        function GroupName = get.GroupName(O)
            GroupName = {O.Data.GroupName}';
        end
    end
    
    %Event listener codes
    methods 
        function handleChangedAx(O, ~, ~)
            %Updates the Gx due to a change in Ax
            Tx = O.Ax; %get Ax since the structure of Gx will be same
            Fields = fieldnames(Tx);
            for f = 1:numel(Fields)
                if ~O.isAxValid(Fields{f}) 
                    Tx.(Fields{f}) = [];
                else
                    Tx.(Fields{f}) = get(Tx.(Fields{f}), 'parent');
                end
            end
            if ~isequal(Tx, O.PrivGx) %Change PrivGx ONLY if there is a change to prevent handleChangedGx triggers.
                O.PrivGx = Tx; 
            end
        end
        
        function handleChangedGx(O, ~, ~)
            %Updates the Ax due to a change in Gx
            Tx = O.Gx; %get Gx since the structure of Ax will be same
            Fields = fieldnames(Tx);
            for f = 1:numel(Fields)
                if ~O.isGxValid(Fields{f})
                    Tx.(Fields{f}) = [];
                else
                    Tx.(Fields{f}) = findobj(Tx.(Fields{f}), 'type', 'axes');
                    if isempty(Tx.(Fields{f}))
                        Tx.(Fields{f}) = axes('Parent', O.Gx.(Fields{f}));
                    end
                end
            end
            if ~isequal(Tx, O.PrivAx) %Change PrivGx ONLY if there is a change to prevent handleChangedGx triggers.
                O.PrivAx = Tx;
            end
        end
    end
    
    %Plot codes
    methods
        function applyAllFigures(O, Method, varargin)
            %Applies a function to all figure handles
            if ~isempty(varargin) && ischar(varargin{1}) && isfield(O.PrivAx, varargin{1})
                Fields = varargin(1);
                varargin = varargin(2:end); %Do only specified axes
            else
                Fields = fieldnames(O.Ax); %Do all axes
            end
            Func = str2func(Method);
            for f = 1:numel(Fields)
                if ~O.isAxValid(Fields{f}); continue; end
                try
                    Func(O.PrivAx.(Fields{f}), varargin{:});
                catch ME
                    throw(ME)
                end
            end
        end
        
        function centerFigureOnMonitor(O, varargin)
            %Centers the figure on the monitor
            varargin = cleanCommandLineInput(varargin{:});
            applyAllFigures(O, 'centerFigureOnMonitor', varargin{:});
        end
        
        function labelSubplots(O, varargin)
            %Labels the panels of subplots
            varargin = cleanCommandLineInput(varargin{:});
            applyAllFigures(O, 'labelSubplots', varargin{:});
        end
        
        function invertFigColor(O, varargin)
            %Sets properties for all axes
            varargin = cleanCommandLineInput(varargin{:});
            applyAllFigures(O, 'invertFigColor', varargin{:});
        end
        
        function setAxes(O, varargin)
            %Sets properties for all axes
            varargin = cleanCommandLineInput(varargin{:});
            if numel(varargin) == 0; return; end
            applyAllFigures(O, 'setAxes', varargin{:});
        end
        
        function setPlotTickDecimal(O, varargin)
            %Sets the number of decimal in the X and/or Y axis
            varargin = cleanCommandLineInput(varargin{:});
            if numel(varargin) == 0; return; end
            applyAllFigures(O, 'setPlotTickDecimal', varargin{:});
        end
        
        function resizeSubplots(O, varargin)
            %Resizes all subplots in every figure
            varargin = cleanCommandLineInput(varargin{:});
            if numel(varargin) == 0; return; end
            applyAllFigures(O, 'resizeSubplots', varargin{:});
        end
        
        function resizeFigure(O, varargin)
            %Resizes all figures
            varargin = cleanCommandLineInput(varargin{:});
            if numel(varargin) == 0
                varargin = {O.Width, O.Height}; %Use the default width and height
            end
            applyAllFigures(O, 'resizeFigure', varargin{:});
        end
        
        function close(O)
            %Closes all figures
            Fields = fieldnames(O.Gx);
            for j = 1:numel(Fields)
                CurGx = O.Gx.(Fields{j});
                if ~isempty(CurGx) && isvalid(CurGx)
                    close(O.Gx.(Fields{j}));
                end
                O.Gx.(Fields{j}) = [];
                O.Ax.(Fields{j}) = [];
            end
        end
        
        function hide(O)
            %Hides all figures that were already made
            Fields = fieldnames(O.Gx);
            for j = 1:numel(Fields)
                CurGx = O.Gx.(Fields{j});
                if ~isempty(CurGx) && isvalid(CurGx)
                    set(O.Gx.(Fields{j}), 'Visible', 'off');
                end
            end
        end
        
        function show(O)
            %Shows all figures that were already made
            Fields = fieldnames(O.Gx);
            for j = 1:numel(Fields)
                CurGx = O.Gx.(Fields{j});
                if ~isempty(CurGx) && isvalid(CurGx)
                    set(O.Gx.(Fields{j}), 'Visible', 'on');
                end
            end
        end
       
        function delete(O)
            %Deletes object and all figures
            close(O)
        end    
        
        function applyDefault(O, varargin)
            %Applies default axes settings to all axes
            applyAllFigures(O, 'setAxes', varargin{:}, O.Default);
%             applyAllFigures(O, 'resizeFigure', varargin{:}, 'FigWidth', O.Width, 'FigHeight', O.Height);
            applyAllFigures(O, 'resizeSubplots', varargin{:});
        end
        
        function setDefault(O, varargin)
            %Sets and applies the default settings to all axes
            for j = 1:2:numel(varargin)
                O.Default.(varargin{j}) = varargin{j+1};
            end
            applyDefault(O);
        end
    end
    
    %Private methods for checking status of an object, like hasData.
    methods (Access = protected)
        function HasData = hasData(O)
            HasData = ~(numel(O.Data) == 0 || (numel(O.Data) == 1 && isempty(O.Data(1).Data)));
        end
    end    
end