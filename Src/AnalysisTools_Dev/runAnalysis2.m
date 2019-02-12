%runAnalysis2 will run the analysis interface for processing single/groups of sequence files.
%
%  runAnlaysis(Analysis, Param, Value)
%
%  Analysis: the following string to perform a specific analysis. If empty, will do all.
%
%  Frequencies
%    VDJfrequency
%    SHMfrequency
%    CDRlength
%    Tree
%
%  Mutations
%
%
%  CDR3 Convergence
%
%  AAScatterPlot
%
%
%



function varargout = runAnalysis2(varargin)
RunInLocalEnv = isempty(varargin);
varargout = cell(1, nargout);

Data = BriliaDatastore(varargin{:});
DataMethods = setdiff(methods('BriliaDatastore'), methods('handle'));
DataFields = fields(Data);

while true
    if RunInLocalEnv
        varargin = input('Analysis> ', 's'); %Start with user input
        if isempty(varargin) 
            continue
        elseif strcmpi(varargin, 'exit')
            return
        end
    end
    varargin = cleanCommandLineInput(varargin);
    
    %Check and correct methods
    try
        %Treat as a method first
        Idx1 = find(strcmpi(varargin{1}, DataMethods), 1);
        if ~isempty(Idx1)
            FH = str2func(DataMethods{Idx1});
            FH(Data, varargin{2:end});
        end
        
        %Treat as a field next
        Idx2 = find(strcmpi(varargin{1}, DataFields), 1);
        if ~isempty(Idx2)
            if numel(varargin) == 2
                Data.(DataFields{Idx2}) = varargin{2};
            else
                disp(Data.(DataFields{Idx2}));
            end
        end
        
        if isempty(Idx1) && isempty(Idx2)
            warning('%s: Unknown methods or fields for BriliaDatastore, "%s".', mfilename, varargin{1});
            fprintf('Valid methods are:\n');
            fprintf('  %s\n', DataMethods{:});
            fprintf('\n');
            fprintf('Valid fields are:\n');
            fprintf('  %s\n', DataFields{:});
        end
    catch ME
        disp(ME)
    end
    if RunInLocalEnv; continue; else; return; end
end