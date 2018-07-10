%joinData will combine a spliced VDJdata from spliceData into a single
%VDJdata. It will ensure each group is assigned a unique group number
%within this data set. 
%  
function VDJdata = joinData(VDJdataS, Map)
if isempty(VDJdataS) || ~iscell(VDJdataS{1})
    error('%s: VDJdataS is not a cell of VDJdata cells.', mfilename);
end

%Renumber each group
VDJdata = vertcat(VDJdataS{:});