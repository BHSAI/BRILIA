function varargout = getIntraDiversity(varargin)
DI = DataFetcher('poll',       'GrpNum, Template', ...
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

VDJdata = spliceData(VDJdata, Map, IdxC);
Indices = zeros(numel(IdxC), 1);
TemplateIdx = Map.Template;
parfor j = 1:numel(IdxC)
    Indices(j) = calcSimpsonDiversity(cell2mat(VDJdata{j}(:, TemplateIdx)));
end
varargout{1} = Indices;