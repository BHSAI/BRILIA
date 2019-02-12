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
function varargout = getHDRF(varargin)
DI = DataFetcher('poll',       'Chain, hGeneName, hLength, hDel, hCDR3', ...
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

%Extract the D gene reading frame
if ~contains(Map.Chain, 'H'); return; end
IdxF = cellfun(@(x) x(1), IdxC);
    
LenVN = sum(cell2mat(VDJdata(IdxF, Map.hLength(1:2))), 2);
DelD5 = cell2mat(VDJdata(IdxF, Map.hDel(2)));
CDR3S = cell2mat(VDJdata(IdxF, Map.hCDR3(3)));
Chart = [1  3  2;   %Each row is the mod(DDel5,3)+1
         2  1  3;   %Each col is the mod(LenVN-CDR3S+1,3)+1
         3  2  1];
R = mod(DelD5, 3) + 1;
C = mod(LenVN-CDR3S+1, 3) + 1;
RF = Chart(sub2ind([3 3], R, C));
DRF = getFirst(VDJdata(IdxF, Map.hGeneName(2)));
UnqD = sort2(unique(DRF));
UnqD = repelem(UnqD, repelem(3, 1, numel(UnqD)));
TRF = 1;
for j = 1:numel(UnqD)
    UnqD{j} = sprintf('%s-RF%d', UnqD{j}, TRF);
    TRF = TRF + 1;
    if TRF > 3; TRF = 1; end
end

parfor j = 1:length(DRF)
    DRF{j} = sprintf('%s-RF%d', DRF{j}, RF(j));
end

%This must be done to ensure RF1-3 are present for all unique D genes. Otherwise, you'll have D1-RF1, D1-RF3, etc.
CurData = countData(DRF);
TempData = countData(UnqD);
ConfData = conformDist(TempData, CurData);
ConfData = ConfData(:, [1 3]);

varargout{1} = ConfData;
