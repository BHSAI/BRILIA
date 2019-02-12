%getTemplate description
%
%  Out = getTemplate(I1, I2, Parameter, Value, ...)
%
%  INPUT
%    I1:
%    I2:
%
%    Parameter       Value      Details
%    --------------- ---------- -------------------------------------------
%
%  OUTPUT
%    TreeTrunk:  matrix of germline to 1st observed ancestor SHM count (not percentage)
%    TreeHeight: matrix of 1st observed ancestor to most mutated descendant SHM count (not percentage)
%
function varargout = getTreeThickness(varargin)
DI = DataFetcher('poll',       'GrpNum, Template, hVmut, hMmut, hDmut, hNmut, hJmut, lVmut, lNmut, lJmut', ...
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

N = numel(IdxC);
TreeThickness  = zeros(N, 1);
MutIdx = nonzeros([Map.hVmut; Map.hMmut; Map.hDmut; Map.hNmut; Map.hJmut; Map.lVmut; Map.lNmut; Map.lJmut]);
for j = 1:N
    TotalSHM = sum(cell2mat(VDJdata(IdxC{j}, MutIdx)), 2);
    TreeThickness(j) = max(TotalSHM) - TotalSHM(1);
end
varargout{1} = TreeThickness;
