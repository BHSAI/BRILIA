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
%    O1:
%    O2:
%
function varargout = getLJGene(varargin)
DI = DataFetcher('poll',       'lGeneName', ...
                 'groupsteps', 'rescale', ...
                 'DataType',   'freq', ...
                 'Level',      'Clonotype');
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

if ~contains(Map.Chain, 'L'); return; end

IdxF = cellfun(@(x) x(1), IdxC);
varargout{1} = getFirst(VDJdata(IdxF, Map.lGeneName(2)));
