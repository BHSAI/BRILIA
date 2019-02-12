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
function varargout = getTreeTwig(varargin)
DI = DataFetcher('poll',       'GrpNum, Template, hRefSeq, hSeq, lRefSeq, lSeq', ...
                 'groupsteps', 'normalize', ...
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

%Extract Twig, parent to child relations
P2C = zeros(size(VDJdata, 1), 1);
IdxDel = cellfun(@(x) x(1), IdxC);
IdxP2C = cellfun(@(x) x(2:end), IdxC, 'un', 0);
IdxP2C = IdxP2C(~cellfun('isempty', IdxP2C));
IdxP2C = vertcat(IdxP2C{:});
for j = 1:numel(IdxP2C)
    TotalSHM = sum(~cmprSeqMEX(VDJdata{IdxP2C(j), Map.hSeq}, VDJdata{IdxP2C(j), Map.hRefSeq}));
    P2C(j) = TotalSHM;
end
P2C(IdxDel) = [];
P2C(P2C == 0) = [];
varargout{1} = P2C;