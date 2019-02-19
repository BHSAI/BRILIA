%DataFetcher serves to extract and process data for many getDataMethods.
%DataFetcher can be used either as a temporary storage object by each
%getData method, or as an object storing static methods for extracting the
%requested data from the BRILIA datastore.
%
%  ReqMapField = getData('poll')
%  GroupSteps = getData('groupsteps')
%  DataType = getData('datatype')
%  Level = getData('level')
%
classdef DataFetcher < handle
    properties 
        Poll        %Fields in Map struct required to collect data
        GroupSteps  %Function to use in Data object to modify data properly
        DataType    %Freq, Pairwise, Lineage, etc..
        Level       %Clonal, Clonotype, CDR3
    end
    
    methods
        function O = DataFetcher(varargin)
            P = inputParser;
            P.KeepUnmatched = 1;
            addParameter(P, 'Poll', {}, @(x) iscell(x) || ischar(x));
            addParameter(P, 'GroupSteps', {}, @(x) iscell(x) || ischar(x));
            addParameter(P, 'DataType', 'Freq', @ischar);
            addParameter(P, 'Level', 'Clonotype', @ischar);
            parse(P, varargin{:})
            O.Poll = P.Results.Poll;
            O.GroupSteps = P.Results.GroupSteps;
            O.DataType = P.Results.DataType;
            O.Level = P.Results.Level;
        end
        
        function set.DataType(O, Type)
            if ischar(Type)
                Type = strsplit(strrep(Type, ' ', ''), ',');
            end
            ValidTypes = {'Freq', 'Corr', 'Tree', 'Confetti', 'Diversity', 'Convergence', 'Lineage', 'Stacked'};
            [~, ValidIdx] = intersect(lower(ValidTypes), lower(Type));
            if isempty(ValidIdx)
                Msg = cell(1, 4);
                Msg{1} = sprintf('%s: Not a valid data type, "%s".', mfilename, Type{1});
                Msg{2} = sprintf(' Valid types are:\n');
                Msg{3} = sprintf(' %s,', ValidTypes{1:end-1});
                Msg{4} = sprintf(' %s.\n', ValidTypes{end});
                warning([Msg{:}]);
            else
                O.DataType = ValidTypes(ValidIdx);
            end
        end
               
        function set.GroupSteps(O, GroupSteps)
            if ischar(GroupSteps)
                O.GroupSteps = strsplit(strrep(GroupSteps, ' ', ''), ',');
            end
            O.GroupSteps = O.GroupSteps(~cellfun('isempty', O.GroupSteps));
        end
        
        function set.Poll(O, Poll)
            Map = getVDJmapper({});
            ValidFields = fieldnames(Map);
            if ischar(Poll)
                Poll = strsplit(strrep(Poll, ' ', ''), ',');
            end
            Poll = Poll(~cellfun('isempty', Poll));
            [~, ValidIdx, PollIdx] = intersect(lower(ValidFields), lower(Poll));
            if numel(PollIdx) ~= numel(Poll) %Some entries were wrong
                PollLoc = ones(size(Poll), 'logical');
                PollLoc(PollIdx) = 0;
                BadFields = Poll(PollLoc);
                Msg = cell(1, 4);
                Msg{1} = sprintf('%s: Some fields are not valid Map fields. These were NOT added:\n', mfilename);
                Msg{2} = sprintf(' %s,', BadFields{1:end-1});
                Msg{3} = sprintf(' %s.\n', BadFields{end});
                warning([Msg{:}]);
            else
                O.Poll = ValidFields(ValidIdx);
            end
        end
        
        function set.Level(O, Level)
            ValidFields = {'', 'Clone', 'Clonotype', 'CDR3', 'AC', 'BC', 'SC'}; %AC is same as Clonotype
            ValidLoc = strcmpi(ValidFields, Level);
            if any(ValidLoc)
                O.Level = ValidFields{find(ValidLoc, 1)};
            else
                Msg = cell(1, 3);
                Msg{1} = sprintf('%s: Specified Level, "%s", is not valid. Valid Levels are:\n', mfilename, Level);
                Msg{2} = sprintf(' %s,', ValidFields{1:end-1});
                Msg{3} = sprintf(' %s.\n', ValidFields{end});
                warning([Msg{:}]);
            end
        end
        
        function [TF, Info] = isAskingForInfo(O, varargin)
            TF = 0;
            Info = [];
            if ~isempty(varargin) && ischar(varargin{1})
                if strcmpi(varargin{1}, 'poll') %this return what Map.(field) is needed to process this data
                    TF = true;
                    Info = O.Poll;
                elseif strcmpi(varargin{1}, 'groupsteps') %this returns what data processing steps are needed to do group comparisons
                    TF = true;
                    Info = O.GroupSteps;
                elseif strcmpi(varargin{1}, 'level')
                    TF = true;
                    Info = O.Level;
                elseif strcmpi(varargin{1}, 'datatype')
                    TF = true;
                    Info = O.DataType;
                end
            end
        end
    end
    
    methods (Static)
        varargout = getHDRF(varargin)
        varargout = getHDGene(varargin)
        varargout = getHJGene(varargin)
        varargout = getHVGene(varargin)
        varargout = getLJGene(varargin)
        varargout = getLVGene(varargin)
        varargout = getMembers(varargin)
        varargout = getTemplate(varargin)
        varargout = getTemplateClonotype(varargin)
        varargout = getTemplateCDR3(varargin)
        varargout = getTreeNodes(varargin)
        varargout = getTreeLeaves(varargin)
        varargout = getTreeTrunk(varargin)
        varargout = getTreeTwig(varargin)
        varargout = getTreeHeight(varargin)
        varargout = getTreeThickness(varargin)
        varargout = getUniqueCDR3(varargin)
        varargout = getIntraDiversity(varargin)
        varargout = getGermCDR3(varargin)
        varargout = getAllCDR3(varargin)
        varargout = getUnqCDR3(varargin)
        varargout = getClonotypeLevel(varargin)
        
        %List the methods
        function Methods = listMethods(varargin)
            %Get all methods available
            AllMethods = setdiff(methods('DataFetcher'), methods('handle'));
            AllMethods = AllMethods(startsWith(AllMethods, 'get'));
            %Return the correct, case-sensitive methods name
            if isempty(varargin)
                Methods = AllMethods;
            else
                if ischar(varargin{1}) && numel(varargin) == 1
                    varargin = strrep(strsplit(varargin{1}, ','), ' ', '');
                end
                %Get the query methods being requested, add those without the leading 'get'
                NoGetStrLoc = ~startsWith(varargin, 'get');
                varargin(NoGetStrLoc) = cellfun(@(x) ['get' x], varargin(NoGetStrLoc), 'un', 0);
                [~, Idx, Idx2] = intersect(lower(AllMethods), lower(varargin)); %case-insensitive serach
                [~, SortIdx] = sort(Idx2);
                Idx = Idx(SortIdx); %To ensure it's in the order of request
                Methods = AllMethods(Idx);
                if isempty(Methods)
                    warning('%s: Could not find method "%s".', mfilename, varargin{1});
                    return
                end
            end
        end
        
        %List the unique fields that is required by all methods
        function Fields = listFields(varargin)
            if nargin == 0
                Methods = DataFetcher.listMethods;
            elseif ischar(varargin{1}) && numel(varargin) == 1
                Methods = strrep(strsplit(varargin{1}, ','), ' ', '');
                Methods = DataFetcher.listMethods(Methods{:});
            else
                Methods = DataFetcher.listMethods(varargin{:});
            end
            Fields = cell(1, numel(Methods));
            for j = 1:numel(Methods)
                Fields{j} = DataFetcher.(Methods{j})('poll');
            end
            Fields = unique(vertcat('Chain', Fields{:})); %Chain must always be there.
        end
        
        %Searches through Map to find the relevant fields being sought, returning a matrix of non-zero indices only
        function MapIdx = listMapIdx(Map, Fields)
            MapIdx = cellfun(@(x) Map.(x), Fields, 'un', 0);
            NumLoc = cellfun(@isnumeric, MapIdx);
            MapIdx = nonzeros(vertcat(MapIdx{NumLoc}));  %isnumeric is used to NOT error out for the Map.Chain = char.
        end
        
        %Format Datastore such that O.SelectedVariableNames includes the minimum number of columns to read to get data
        function prepDatastoreForFetch(O, varargin)
            if nargin == 1
                GetMethods = DataFetcher.listMethods;
            elseif ischar(varargin{1}) && numel(varargin) == 1
                GetMethods = strrep(strsplit(varargin{1}, ','), ' ', '');
            else
                GetMethods = varargin;
            end
            GetFields  = DataFetcher.listFields(GetMethods{:});
            GetMapIdx  = DataFetcher.listMapIdx(O.Map, GetFields);
            O.SelectedVariableNames = GetMapIdx;
        end
        
        %Let the DataFetcher fetch the data that is requested from all files AND assemble them into the object
        function fetchData(O, varargin)
            %Determine the Level - AC, BC, SC, etc
            LevelIdx = find(strcmpi(varargin, 'level'), 1);
            if ~isempty(LevelIdx)
                Level = upper(varargin{LevelIdx+1});
                varargin(LevelIdx:LevelIdx+1) = [];
            else
                Level = 'AC';
            end
            switch Level
                case {'AC', 'CLONOTYPE'}
                    Idx = O.Idx;
                case 'BC'
                    Idx = O.Idx;
                    for j = 1:numel(O.Idx)
                        Idx{j} = Idx{j}(cellfun('length', Idx{j}) >  1);
                    end
                case 'SC'
                    Idx = O.Idx;
                    for j = 1:numel(O.Idx)
                        Idx{j} = Idx{j}(cellfun('length', Idx{j}) == 1);
                    end
                case 'CLONE'
                    Idx = O.Idx;
                    for j = 1:numel(O.Idx)
                        Idx{j} = num2cell(vertcat(Idx{j}{:}));
                    end
            end
            
            GetMethods = DataFetcher.listMethods(varargin{:});       %Get the specified, or all, data get methods
            if isempty(GetMethods); return; end
            DataNames = cellfun(@(x) x(4:end), GetMethods, 'un', 0); %DataNames are the methods without the leading "get"
            DataFetcher.prepDatastoreForFetch(O, GetMethods{:}); %Ensures the SelectedVariableNames field is set correctly
            
            %Fetch the raw data per file
            FetchedData(1:numel(O.Files)) = struct;
            O.reset;
            for f = 1:numel(O.Files)
                fprintf('%s: File %d of %d.\n', mfilename, f, numel(O.Files));
                [VDJdata, Map] = read(O);
                for m = 1:numel(GetMethods)
                    FetchedData(f).(DataNames{m}) = DataFetcher.(GetMethods{m})(VDJdata, Map, Idx{f});
                end
            end
            
            %Place each data into the right DataObject
            Fields = fieldnames(FetchedData);
            DataObj = cell(1, numel(Fields));
            q = 1;%Object's qth number
            for f = 1:numel(Fields)
                if all(cellfun('isempty', {FetchedData.(Fields{f})}))
                    continue
                end
                DataType = DataFetcher.(GetMethods{f})('DataType');
                GroupSteps = DataFetcher.(GetMethods{f})('GroupSteps');
                for d = 1:numel(DataType)
                    switch lower(DataType{d})
                        case 'freq' %plots bar graphs
                            ExtData = {FetchedData.(Fields{f})};
                            if isstruct(ExtData{1}) %This one outputs a structure with field S.Data, S.Idx, etc. take out S.Data.
                                ExtData = cellfun(@(x) x.Data, ExtData, 'un', 0);
                            end
                            DataObj{q} = FrequencyData(FetchedData.(Fields{f}), 'Group', O.Group, 'DataName', Fields{f}, 'CtrlGroup', 1); %CtrlGroup is always 1
                            for g = 1:numel(GroupSteps)
                                Steps = strsplit(GroupSteps{g}, {'=', '-', ':'});
                                Steps = cleanCommandLineInput(Steps{:});
                                DataObj{q}.(Steps{1})(Steps{2:end});
                            end
                            q = q + 1;
%                         case 'diversity' %Plots confetti map, entropy bar, diversity bar, template scatters
%                             DataObj{q} = DiversityData({FetchedData.(Fields{f})}, 'Group', O.Group, 'DataName', Fields{f}, 'CtrlGroup', 1); %CtrlGroup is always 1
%                             q = q + 1;
                           case 'convergence' %Plots NxN correlation matrix for CDR3, bar of convergent CDR3s, dendrogram, etc
                              DataObj{q} = ConvergenceData(FetchedData.(Fields{f}), 'Group', O.Group, 'DataName', Fields{f}, 'CtrlGroup', 1); %CtrlGroup is always 1
                              q = q + 1;
%                         case 'stacked'
%                             DataObj{q} = StackedData({FetchedData.(Fields{f})}, 'Group', O.Group, 'DataName', Fields{f}, 'CtrlGroup', 1); %CtrlGroup is always 1
%                             q = q + 1;
                    end
                end
            end
            O.set('DataObj', DataObj(~cellfun('isempty', DataObj))); %Removes failed/empty DataObj
        end
    end
end