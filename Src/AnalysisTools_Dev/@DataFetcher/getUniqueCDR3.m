%getTemplate description
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
function varargout = getUniqueCDR3(varargin)
% CODING_NOTE: Modify this for the data being collected
DI = DataFetcher('Poll',       'Chain, GrpNum, hCDR3, lCDR3', ...
                 'GroupSteps', 'rescale', ...
                 'DataType',   'freq', ...
                 'Level',      'clonotype');

% CODING_NOTE: Do NOT modify this input parsing
[TF, varargout{1}] = DI.isAskingForInfo(varargin{:});
if TF; return; end
if nargin == 2
    [VDJdata, Map] = deal(varargin{1:2});
    IdxC = getGrpIdx(VDJdata, Map, 'AC');
elseif nargin == 3
    [VDJdata, Map, IdxC] = deal(varargin{1:3});
else
    error('%s: Incorrect number of inputs.', mfilename);
end

% CODING_NOTE: Modify below to for data collection
VDJdata = spliceData(VDJdata, Map, IdxC);
N = numel(IdxC);
Out = zeros(N, numel(Map.Chain));
Idx = nonzeros([Map.hCDR3(1); Map.lCDR3(1)]);
parfor j = 1:N
    CDR3 = VDJdata{j}(:, Idx);
    Out(j) = numel(unique2(CDR3), 'rows');
end
varargout{1} = Out;
