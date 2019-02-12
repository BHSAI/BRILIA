%getAllCDR3 description
%
%  [O1, O2] = getTemplate(I1, I2, Parameter, Value, ...)
%
%  INPUT
%    I1:
%    I2:
%
%    Parameter       Value      Details
%    --------------- ---------- -------------------------------------------
%
%  OUTPUT
%    Out: frequency of unique hCDR3/lCDR3/hCDR3+lCDR3 per clonotype
%
function varargout = getAllCDR3(varargin)
% CODING_NOTE: Modify this for the data being collected
DI = DataFetcher('Poll',       'Chain, GrpNum, Template, hCDR3, hRefSeq, lCDR3, lRefSeq', ...
                 'GroupSteps', '', ...
                 'DataType',   'convergence', ...
                 'Level',      'clone');

% CODING_NOTE: Do NOT modify this input parsing
[TF, varargout{1}] = DI.isAskingForInfo(varargin{:});
if TF; return; end
if nargin >= 2
    [VDJdata, Map] = deal(varargin{1:2});
else
    error('%s: Incorrect number of inputs.', mfilename);
end

% CODING_NOTE: Modify below to for data collection
for c = 1:numel(Map.Chain)
    C = lower(Map.Chain(c));
    Out.([C 'CDR3']) = VDJdata(:, Map.([C 'CDR3'])(1));
end
Out.GrpNum = cell2mat(VDJdata(:, Map.GrpNum));
Out.Template = cell2mat(VDJdata(:, Map.Template));

varargout{1} = Out;
