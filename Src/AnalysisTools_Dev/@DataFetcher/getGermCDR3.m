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
function varargout = getGermCDR3(varargin)
% CODING_NOTE: Modify this for the data being collected
DI = DataFetcher('Poll',       'Chain, GrpNum, Template, hCDR3, hRefSeq, lCDR3, lRefSeq', ...
                 'GroupSteps', '', ...
                 'DataType',   'convergence', ...
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
GrpIdx = cellfun(@(x) x(1), IdxC);
for c = 1:numel(Map.Chain)
    C = lower(Map.Chain(c));
    %Get the RefSeq, then determine where the CDR3 is before converting
    Out.([C 'CDR3']) = cell(numel(GrpIdx), 1);
    for j = 1:numel(GrpIdx)
        RefSeq = VDJdata{GrpIdx(j), Map.([C 'RefSeq'])};
        CDR3S  = VDJdata{GrpIdx(j), Map.([C 'CDR3'])(3)};
        CDR3E  = VDJdata{GrpIdx(j), Map.([C 'CDR3'])(4)};
        Out.([C 'CDR3']){j} = RefSeq(CDR3S:CDR3E); %Store nt sequences and translate to AA all at once below.
    end
    Out.([C 'CDR3']) = nt2aa(Out.([C 'CDR3']), 'acgtonly', false);
end
Out.GrpNum = cell2mat(VDJdata(GrpIdx, Map.GrpNum));
Out.Template = regroupTemplate(cell2mat(VDJdata(:, Map.Template)), cell2mat(VDJdata(:, Map.GrpNum)));

varargout{1} = Out;
