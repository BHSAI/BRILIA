%getClonotypeLevel description
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
function varargout = getClonotypeLevel(varargin)
DI = DataFetcher('poll',       'GrpNum, Template', ...
                 'groupsteps', 'normalize', ...
                 'DataType',   'stacked', ... 
                 'Level',      'Clonotype');
[TF, varargout{1}] = DI.isAskingForInfo(varargin{:});
if TF; return; end

if nargin == 2
    [VDJdata, Map] = deal(varargin{1:2});
    IdxC = getGrpIdx(VDJdata, Map, 'AC');
elseif nargin == 3
    [VDJdata, Map, ~] = deal(varargin{1:3});
    [~, ~, ~, IdxC] = unique2(VDJdata(:, Map.GrpNum));
else
    error('%s: Incorrect number of inputs.', mfilename);
end



%Extract AC,BC,SC numbers
AC = numel(IdxC);
BC = sum(cellfun('length', IdxC) > 1);
SC = AC-BC;
varargout{1} = {'SC' SC;...
                'BC' BC};

% Member = cellfun('length', IdxC);
% Template = zeros(numel(IdxC), 1);
% for j = 1:numel(Template)
%     Template(j) = sum(cell2mat(VDJdata(IdxC{j}, Map.Template)));
% end% 
% %For a cutoff of TC <= 5 for Small clonotypes
% Cutoff = 5;
% BC_Small = sum(Member > 1 & Template <= Cutoff);
% BC_Large = sum(Member > 1 & Template >  Cutoff); 
% SC_Small = sum(Member == 1 & Template <= Cutoff);
% SC_Large = sum(Member == 1 & Template >  Cutoff);
% 
% varargout{1} = {'BC_Small' BC_Small;...
%                 'BC_Large' BC_Large;...
%                 'SC_Small' SC_Small;...
%                 'SC_Large' SC_Large};
