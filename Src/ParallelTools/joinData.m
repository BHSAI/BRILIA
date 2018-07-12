%joinData will combine a spliced VDJdata from spliceData into a single
%VDJdata. It WILL renumber groups to ensure each group number is unique. To
%join a spliced VDJdata without renumbering groups, simply use vertcat. 
%
%  INPUT
%    VDJDataS: spliced VDJdata that stores a cell with cell arrays of
%      VDJdata split by group number
%    Map: strcuture of indices of VDJdata cell headers
%   
%  OUTPUT
%    VDJdata: a MxN cell array of VDJdata
%
%  EXAMPLE
%    [VDJdata, VDJheader] = openSeqData(); %Choose a BRILIA output file
%    VDJdata = spliceData(VDJdata, Map)
%    parfor j = 1:length(VDJdata)
%       %DO SOMETHING BY GROUP
%    end
%    VDJdata = joinData(VDJdata); 
%  
function VDJdataS = joinData(VDJdataS, Map)
if isempty(VDJdataS) || ~iscell(VDJdataS{1})
    error('%s: VDJdataS is not a cell of VDJdata cells.', mfilename);
end
if iscell(Map)
    Map = getVDJmapper(Map);
end

GrpNumShift = 0;
for j = 1:length(VDJdataS)
    if isempty(VDJdataS{j}); continue; end
    VDJdataS{j} = renumberVDJdata(VDJdataS{j}, Map, 'grp');
    GrpNum = cell2mat(VDJdataS{j}(:, Map.GrpNum)) + GrpNumShift;
    GrpNumShift = max(GrpNum);
    VDJdataS{j}(:, Map.GrpNum) = num2cell(GrpNum);
end
VDJdataS = vertcat(VDJdataS{:});