%getProdVDJ will take only the productive VDJ's from the VDJdata file. A
%productive VDJ data follows these rules. Based on results provided by
%LabelNonprodVDJ.m, where "N" is nonproductive. "M and Y" are accepted as
%productive here.
%  1) No stop codon
%  2) No pseudogene for V
%
%  VDJdata = getprodVDJ()
%  VDJdata = getprodVDJ(VDJdata,NewHeader)
%  [VDJdata, VDJdataNP] = getprodVDJ(...)    will return the nonprodVDJ's
%  too.

function [VDJdata,varargout] = getProdVDJ(varargin)
if isempty(varargin)
    [VDJdata, NewHeader, ~, ~] = openSeqData;
else
    VDJdata = varargin{1};
    NewHeader = varargin{2};
end
getHeaderVar;

KeepThis = ones(size(VDJdata,1),1,'logical');
for j = 1:length(size(VDJdata,1))
    if strcmpi(VDJdata{j,FunctLoc},'N')
        KeepThis(j) = 0;
    end
end

if nargout == 2
    varargout{1} = VDJdata(KeepThis == 0,:);
end
VDJdata = VDJdata(KeepThis,:);