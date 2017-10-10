%joinData will combine a split VDJdata from splitDataByGroup into a single
%VDJdata. It will ensure each group is assigned a unique group number
%within this data set. If combining with another ParVDJdata, must update
%the group numbers.

function VDJdata = joinData(ParVDJdata)
if isempty(ParVDJdata) || ~iscell(ParVDJdata{1})
    error('%s: ParVDJdata is not a cell of VDJdata cells.', mfilename);
end
VDJdata = cat(1, ParVDJdata{:});