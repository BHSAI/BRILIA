%Datastore is used to track which BRILIA output files to use in the
%analysis, and which files are grouped according to the same treatment.
%
%  ATTRIBUTES
%    DS          Datastore for accessing data directly from hard drive
%    Map         Structure mapping a name to column number
%    VDJheader   Cell array of unmodified header names
%    Idx         Cell array of cell arrays of all clonotype indices
%    Files       Cell array of BRILIA file names
%    Group       Int  array of treatment groups of size(Files)
%    GroupName   Cell array of unique group name of size(unique(Group))
%
%  CONSTRUCTOR
%    Datastore(FileNames, 'GroupBy', GroupBy, 'GroupFormat', 'GroupFormat)
%      FileNames: cell array of file names. If empty, will ask user.
%      GroupBy ['file' 'dir']: Group files by treatment group based on a 
%        unique file name or folder name.
%      GroupFormat: regular expression to search to group file or folders. 
%        Ex: 'grp#' will group data by file/folder name containing the 
%        "grp1" together, "grp2" together, etc.
%
%  METHODS
%    [Data, Map, VariableNames] = read(O, SelectedVariableNames)
%      Reads the datastore at SelectedVariableNames (by name or column
%      number), returning cell array of data, Map structure, and
%      VariableNames.
%
%    [Data, Map, VariableNames] = readall(O, SelectedVariableNames)
%      Reads the datastore at SelectedVariableNames (by name or column
%      number), returning cell array of cell arrays of data, Map structure,
%      and VariableNames.
%
%    groupBy(O, Options, Format)
%      Sets the file grouping by Options ['file', 'dir'] and using regular
%      expression defined by Format, such as 'Grp#'.
%
%  NOTE
%    When handling multiple data files, ALL data files MUST have the SAME
%    number of columns, in the SAME order. This is not an issue for
%    individual files.

classdef Datastore < handle
    properties (AbortSet, Access = public, Dependent)
        SelectedVariableNames %This is used by the datastore
        Group            %Group number
        Files            %Cell array of BRILIA file names
    end
    
    properties (AbortSet, GetAccess = public, SetAccess = private)
        DS                %Datastore for accessing data directly from hard drive
        Map               %Structure mapping a name to column number
        VDJheader         %Cell array of unmodified header names
        Idx               %Cell array of cell arrays of all clonotype indices
        GroupName         %Cell array of unique group name of size(unique(Group))
        CtrlGroupName     %Control group name or number
        Diversity struct  %Graphics array of clonotype data of all files
        DataObj cell      %Struct of cell array data of frequency data
    end
    
    properties (Access = public)
        SaveDir char     %Parent directory to save all outputs
    end
    
    methods
        function O = Datastore(FileNames, varargin)
            P = inputParser;
            addParameter(P, 'GroupBy',      'file', @(x) ischar(x) && ismember(lower(x), {'file', 'dir', 'folder'}));
            addParameter(P, 'GroupFormat',  'Grp#', @ischar);
            parse(P, varargin{:});
            GroupBy = P.Results.GroupBy;
            GroupFormat = P.Results.GroupFormat;
            
            if nargin < 1 || isempty(FileNames)
                TmpFiles = getBriliaFiles('', true, false);
                assert(~isempty(TmpFiles), '%s: No BRILIA files were chosen for the datastore.', mfilename);
            elseif isa(FileNames, 'Datastore')
                TmpFiles = FileNames.Files;
            else
                TmpFiles = FileNames;
            end
            
            %Setup the datastore properties
            O.DS = datastore(TmpFiles);
            O.DS.ReadSize = 'file';
            O.VDJheader = readDlmFile(O.DS.Files{1}, 'LineRange', 1);
            O.Map = getVDJmapper(O.VDJheader);
            O.fixTextscanFormats; %Needs to be prevent error reading char as a number.
            O.Idx = cell(length(O.Files), 1); %Get the group indices once
            for f = 1:numel(O.Files)
                GrpNum = read(O, O.Map.GrpNum); 
                [~, ~, ~, O.Idx{f}] = unique2(GrpNum);
            end
            O.SelectedVariableNames = ''; %Reset to default
            O.reset; %Reset to beginning
            O.groupBy(GroupBy, GroupFormat); %Assign group name and number
        end
    
        %Read the data file-by-file, starting after last read file. 
        function [Data, Map, VariableNames] = read(O, varargin)
            [Data, Map, VariableNames] = deal([]);
            if ~O.DS.hasdata
                fprintf('%s: End of file reached. Use "reset" to go to the beginning.\n', mfilename);
                return
            end
            if ~isempty(varargin)
                O.SelectedVariableNames = varargin;
            end
            Data = table2cell(read(O.DS));
            if nargout >= 2
                VariableNames = O.SelectedVariableNames;
                Map = getVDJmapper(VariableNames);
                Map.Chain = O.Map.Chain; %Required as Map.Chain is special, and a char.
            end
        end
        
        %Read all data, starting from the beggining
        function [Data, Map, VariableNames] = readall(O, varargin)
            O.SelectedVariableNames(varargin{:});
            reset(O.DS);
            Data = cell(size(O.DS.Files));
            for j = 1:size(Data, 1)
                Data{j} = table2cell(read(O.DS));
            end
            reset(O.DS);
            if nargout >= 2
                VariableNames = O.SelectedVariableNames;
                Map = getVDJmapper(VariableNames);
                Map.Chain = O.Map.Chain; %%Required as Map.Chain is special, and a char.
            end
        end
        
        %Go back to the beginning
        function O = reset(O)
            reset(O.DS);
        end
        
        %Joins multiple data files into a single file
        O = joinFiles(O, Option, Level);
    end
    
    %Get and set methods
    methods
        function PropVal = get.SelectedVariableNames(O)
            PropVal = O.DS.SelectedVariableNames;
        end
        
        function set.SelectedVariableNames(O, varargin)
            if isempty(varargin{1}) %set to default
                O.DS.SelectedVariableNames = O.DS.VariableNames;
            else
                if iscell(varargin{1}) %Unwrap cell
                    varargin = varargin{1};
                end
                if all(cellfun('isclass', varargin, 'char')) %Directly use name
                    O.DS.SelectedVariableNames = varargin;
                elseif all(cellfun('isclass', varargin, 'double')) %Need to convert number to name
                    O.DS.SelectedVariableNames = O.DS.VariableNames(cell2mat(varargin));
                else
                    error('%s: All inputs must be of same type, char or double.', mfilename);
                end
            end
        end
        
        function set.Files(O, Files)
            O.DS.Files = Files;
        end
        
        function Files = get.Files(O)
            Files = O.DS.Files;
        end
        
        function set.SaveDir(O, SaveDir)
            if nargin < 2 || isempty(SaveDir)
                SaveDir = uiputdir2('*.png', 'Select dir to save plots');
            end
            O.SaveDir = SaveDir;
            if ~isdir(O.SaveDir)
                fprintf('%s: "%s" does does not exist yet.\nWill attempt to create it when needed.\n', mfilename, O.SaveDir);
            end
        end
        
        function setSaveDir(O, SaveDir)
            O.SaveDir = SaveDir;
        end
        
        function O = set(O, Param, Value)
            O.(Param) = Value;
        end
        
        function Value = get(O, Param)
            Value = O.(Param);
        end
    end
    
    %Group assignment methods
    methods 
        function set.Group(O, Group)
            if nargin < 2
                O.GroupName = 1:numel(O.Files);
            else
                if numel(Group) == 1
                    O.GroupName = repelem(Group, 1, numel(O.Files));
                elseif numel(Group) == numel(O.Files)
                    O.GroupName = Group;
                else
                    error('%s: Group must be equal to number of files.', mfilename);
                end
            end
            O.sortGroup;
        end
        
        function Group = get.Group(O)
            [~, ~, Group] = unique(O.GroupName, 'stable');
        end
        
        function O = groupBy(O, Options, Format)
            %Determine the grouping name and number
            if strcmpi(Options, 'file') %GroupName by file name
                GroupNameTmp = O.Files;
                for j = 1:numel(GroupNameTmp)
                    [~, ~, ~, GroupNameTmp{j}] = parseFileName(GroupNameTmp{j});
                end
            elseif isnumeric(Options) %Group by number (useful for random group assignment)
                if numel(Options) == 1
                    TmpGroup = repelem(Options, 1, numel(O.DS.Files));
                    O.GroupName = cellfun(@num2str, (num2cell(TmpGroup)), 'un', 0);
                elseif numel(Options) == numel(O.DS.Files)
                    TmpGroup = O.Options;
                    O.GroupName = cellfun(@num2str, (num2cell(TmpGroup)), 'un', 0);
                else
                    error('%s: The number of group assignments (%d) and files (%d) are unequal.', mfilename, numel(O.Options), numel(O.DS.Files));
                end
                return
            else %GroupName by folder name, 2-level above
                GroupNameTmp = O.Files;
                for j = 1:numel(GroupNameTmp)
                    SlashIdx = find(GroupNameTmp{j} == filesep, 3, 'last');
                    GroupNameTmp{j} = GroupNameTmp{j}(SlashIdx(1)+1:SlashIdx(2)-1);
                end
            end
            if nargin < 3 || isempty(Format) || ~contains(Format, '#')
                O.GroupName = GroupNameTmp;
                return
            end
            Format = strrep(Format, '#', '([a-zA-Z0-9]+)'); %Replace # with the word search
            if ~isempty(Format)
                for j = 1:numel(GroupNameTmp)
                    GroupChar = regexp(GroupNameTmp{j}, Format, 'tokens');
                    if ~isempty(GroupChar) %Unwrap regexp output
                        GroupNameTmp{j} = GroupChar{1}{1};
                    end
                end
            end
            O.GroupName = GroupNameTmp;
        end
        
        function O = sortGroup(O, SortIdx)
            if nargin < 2
                [O.GroupName, SortIdx] = sort(O.Group);
                if all(diff(SortIdx) == 1); return; end %No change
                O.Files = O.Files(SortIdx);
                O.fixTextscanFormats; %This MUST be done as re-ordering DS will change the textscan formats!
            else
                if all(diff(SortIdx) == 1); return; end %No change
                if numel(O.GroupName) == numel(unique(SortIdx))
                    O.Files = O.Files(SortIdx); 
                    O.fixTextscanFormats; %This MUST be done as re-ordering DS will change the textscan formats!
                    O.Idx = O.Idx(SortIdx);
                    O.GroupName = O.GroupName(SortIdx);
                else
                    warning('%s: to use sortGroup, must input a sort index of 1 to %d.', mfilename, numel(O.Files));
                    return
                end
            end
        end
        
        function O = setCtrlGroup(O, CtrlGroupName)
            if isnumeric(CtrlGroupName) && ~isnumeric(O.GroupName)
                CtrlGroupName = num2str(CtrlGroupName);
                CtrlLoc = strcmpi(CtrlGroupName, O.GroupName);
                if ~any(CtrlLoc)
                    warning('%s: No match to CtrlGroupName "%s"', mfilename, CtrlGroupName);
                    return
                end
            elseif isnumeric(CtrlGroupName) && isnumeric(O.GroupName)
                CtrlLoc = CtrlGroupName == O.GroupName;
                if ~any(CtrlLoc)
                    warning('%s: No match to CtrlGroupName "%f"', mfilename, CtrlGroupName);
                    return
                end
            elseif iscell(O.GroupName) && ischar(CtrlGroupName)
                CtrlLoc = strcmpi(CtrlGroupName, O.GroupName);
                if ~any(CtrlLoc)
                    warning('%s: No match to CtrlGroupName "%s"', mfilename, CtrlGroupName);
                    return
                end
            else
                warning('%s: Unable to match to CtrlGroupName type %s (must be char or number) to GroupName type %s (must be cell or matrix).', mfilename, class(CtrlGroupName), class(O.GroupName));
                return
            end
            [~, ~, UnqIdx] = unique(O.GroupName);
            UnqIdx(CtrlLoc) = -1;
            [~, SortIdx] = sort(UnqIdx);
            O.sortGroup(SortIdx);
            O.CtrlGroupName = CtrlGroupName;
        end
    end
    
    %Methods of data handling
    methods
        function fetchData(O, varargin)
            DataFetcher.fetchData(O, varargin{:});
            if ~isempty(O.SaveDir)
                for j = 1:numel(O.DataObj)
                    O.DataObj{j}.SaveDir = O.SaveDir;
                end
            end
        end
        
        %Clears all data and restarts
        function clearData(O)
            O.DataObj = {};
        end
        
        %Plots all clonotype frequency data
        function plotFreqData(O, varargin)
            for j = 1:numel(O.DataObj)
                O.DataObj{j}.setError('minmax');
                O.DataObj{j}.plot(varargin{:});
            end
        end
        
        %Closes all plots
        function close(O)
            for j = 1:numel(O.DataObj)
                O.DataObj{j}.close;
            end
        end
               
        %Closes all plots before deleting Datastore
        function delete(O)
            O.close;
        end
        
        %Saves the underlying data in SaveDir/DataName.csv
        function saveData(O, Field, SaveDir)
            
            if nargin < 2 || isempty(Field)
                warning('%s: No field was specified, "%s"', mfilename, Field);
                return
            elseif ~any(strcmpi(fieldnames(O), Field))
                warning('%s: Not a valid data field to save, "%s"', mfilename, Field);
                return
            end
            ValidFields = fieldnames(O);
            ValidLoc = strcmpi(ValidFields, Field);
            Field = ValidFields(ValidLoc);
            if isempty(O.SaveDir) && (nargin < 3 || isempty(SaveDir))
                SaveDir = uiputdir2('', 'Select dir to save plots');
                if isempty(SaveDir)
                    warning('%s: No save directory was chosen. Aborting save.', mfilename);
                    return
                end
                O.SaveDir = SaveDir;
            end
            [O, Success] = prepSaveDir(O);
            if ~Success; return; end
            for k = 1:numel(Field)
                if isstruct(O.(Field{k})) %Need to save each subfield
                    DataField = fieldnames(O.(Field{k}));
                    if ~isstruct(O.(Field{k})(1).(DataField{1}))
                        Data = struct2array(O.(Field{k}));
                        writeDlmFile(Data, fullfile(O.SaveDir, sprintf('%s.csv', Field{k})));
                    else
                        for j = 1:numel(DataField)
                            Data = O.(Field{k}).(DataField{j});
                            if isstruct(Data)
                                Data = struct2array(Data);
                            end
                            writeDlmFile(Data, fullfile(O.SaveDir, sprintf('%s.%s.csv', Field{k}, DataField{j})));
                        end
                    end
                else
                    writeDlmFile(O.(Field{k}), fullfile(O.SaveDir, sprintf('%s.csv', Field{k})));
                end
            end
        end
        
        %Save the plots in SaveDir/DataName.png
        function savePlot(O, SaveDir)
            if isempty(O.SaveDir) && (nargin < 2 || isempty(SaveDir))
                SaveDir = uiputdir2('*.png', 'Select dir to save plots');
                if isempty(SaveDir)
                    return
                end
                O.SaveDir = SaveDir;
            end
            [O, Success] = prepSaveDir(O);
            if ~Success; return; end
            
            for j = 1:numel(O.DataObj)
                if ~isempty(O.DataObj{j}.Gx) && isvalid(O.DataObj{j}.Gx)
                    O.DataObj{j}.OutputDir = O.SaveDir;
                    O.DataObj{j}.savePlot;
                end
            end
        end
        
        %Lists which plots are available
        function listData(O)
            ObjName = cellfun(@(x) x.DataName, O.DataObj, 'un', 0);
            ObjType = cellfun(@(x) class(x), O.DataObj, 'un', 0);
            ObjNameType = cellfun(@(x, y) [x ' (' y ')'], ObjName, ObjType, 'un', 0);
            if ~isempty(ObjNameType)
                dispList(ObjNameType)
            else
                fprintf('%s: No DataObj to list.\n', mfilename);
            end
        end
        
        %Returns the DataObj that was specified
        function Obj = getData(O, varargin)
            if nargin == 1
                listData(O);
                Obj = [];
                return
            end
            Obj = O.DataObj(findDataObj(O, varargin{:}));
            if numel(Obj) == 1
                Obj = Obj{1};
            end
        end
        
        %Shows a specific data plot
        function plotData(O, DataNameOrNum)
            if nargin == 1
                for j = 1:numel(O.DataObj)
                    O.DataObj{j}.plot;
                end
            else
                ObjNum = findDataObj(O, DataNameOrNum);
                for j = 1:numel(ObjNum)
                    O.DataObj{ObjNum}.plot;
                end
            end
        end
    end
    
    methods (Access = private)
        %Prepares the directory on hard drive for
        function [O, Success] = prepSaveDir(O)
            Success = true;
            if isempty(O.SaveDir)
                warning('%s: SaveDir is undefined. Use set(OBJ, ''SaveDir'', SAVE_DIR) to set SaveDir', mfilename);
                Success = false;
                return
            end
            if ~isdir(O.SaveDir)
                [Success, ErrMsg] = mkdir(O.SaveDir);
                if ~Success
                    warning('%s: Could not make the output directory, "%s".\nResetting SaveDir to empty.\nError: %s\n', mfilename, O.SaveDir, ErrMsg);
                    O.SaveDir= [];
                end
                return
            end
        end
        
        %Ensures that the MapNum columns are char for datastore reading
        function O = fixTextscanFormats(O)
            CharIdx = nonzeros([O.Map.hGeneNum; O.Map.lGeneNum; O.Map.hOverSeq5; O.Map.hOverSeq3; O.Map.lOverSeq5; O.Map.lOverSeq3]);
            O.DS.TextscanFormats(CharIdx) = {'%q'}; %Correct the GeneNum "[# # #]" formats to Str. This will no longer be used in future releases.
        end
        
        %Find the number of the DataObj based on data name or number
        function ObjNum = findDataObj(O, DataNameOrNum)
            if isnumeric(DataNameOrNum)
                ObjNum = DataNameOrNum;
            else
                ObjNum = find(strcmpi(cellfun(@(x) x.DataName, O.DataObj, 'un', 0), DataNameOrNum));
            end
            if isempty(ObjNum) || min(ObjNum) <= 0 || max(ObjNum) > numel(O.DataObj)
                ObjNum = [];
                fprintf('%s: No valid object was selected. Valid choices are: \n', mfilename);
                listData(O);
            end
        end
    end
end