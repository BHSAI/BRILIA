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
function varargout = getTreeNodes(varargin)
DI = DataFetcher('poll',       'GrpNum, Template, ChildCount', ...
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

%Extract TreeNodes and TreeLeaves
N = numel(IdxC);
TreeNodes = zeros(N, 1);
ChildCount = cell2mat(VDJdata(:, Map.ChildCount));
for j = 1:N
    TreeNodes(j) = sum(ChildCount(IdxC{j}) > 0);
end
varargout{1} = TreeNodes;

