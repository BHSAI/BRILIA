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
function varargout = getTemplateCDR3(varargin)
% CODING_NOTE: Modify this for the data being collected
DI = DataFetcher('Poll',       'Chain, Template, hCDR3, lCDR3', ...
                 'GroupSteps', 'normalize, setError=minmax', ...
                 'DataType',   'freq, diversity', ...
                 'Level',      '');

% CODING_NOTE: Do NOT modify this input parsing
[TF, varargout{1}] = DI.isAskingForInfo(varargin{:});
if TF; return; end
if nargin == 2
    [VDJdata, Map] = deal(varargin{1:2});
elseif nargin == 3
    [VDJdata, Map, ~] = deal(varargin{1:3});
else
    error('%s: Incorrect number of inputs.', mfilename);
end

% CODING_NOTE: Modify below for data collection
Template = cell2mat(VDJdata(:, Map.Template));
CDR3Idx = nonzeros([Map.hCDR3(1); Map.lCDR3(1)]);
[Template, Idx] = regroupTemplate(Template, VDJdata(:, CDR3Idx)); %sum template by unique CDR3 combos
varargout{1} = struct('Data', Template, 'Idx', {Idx}); %Struct form since multiple output required