%spliceData will split VDJdata into smaller cells of VDJdata to enable
%parallel processing of groups.
%
%  VDJdataS = spliceData(VDJdata, Map)
%
%  INPUT
%    VDJdata: BRILIA data cell
%    Map: BRILIA header cell or structure map of header locations
%
%  OUTPUT
%    VDJdataS: Mx1 cell of VDJdata(GrpIdx, :) cell arrays, used to put each
%      group in its own cell for parallel processing.
%
%  EXAMPLE
%    [VDJdata, VDJheader] = openSeqData(); %Choose a BRILIA output file
%    VDJdata = spliceData(VDJdata, Map)
%    parfor j = 1:length(VDJdata)
%       %DO SOMETHING BY GROUP
%    end
%    VDJdata = joinData(VDJdata); 
%    
function VDJdataS = spliceData(VDJdata, Map, IdxC)
if isempty(VDJdata)
    VDJdataS = VDJdata; 
    return 
end
Map = getVDJmapper(Map);
if iscell(VDJdata{1}) %Already spliced
    VDJdataS = VDJdata;
    return
end
if nargin == 2
    [~, ~, ~, IdxC] = unique2(VDJdata(:, Map.GrpNum));
end
VDJdataS = cell(numel(IdxC), 1);
for k = 1:numel(IdxC)
    VDJdataS{k} = VDJdata(IdxC{k}, :);
end