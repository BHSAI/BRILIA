%getVDJmapper returns a structure where the fields refer to column
%number(s) corresponding to a data type in VDJdata cell array. It uses the
%VDJheader information.
%
%  Map = getVDJmappper(VDJheader)
%
%  INPUT
%    VDJheader: 1xN cell of header string for VDJdata
%
%  OUTPUT
%    Map: 

function Map = getVDJmapper(VDJheader)
if isstruct(VDJheader)
    Map = VDJheader;
    return
end
[H, L, Chain] = getAllHeaderVar(VDJheader);
CommonFields = {'SeqNameLoc', 'SeqNumLoc', 'GrpNumLoc', 'TemplateLoc', 'ChildCountLoc'};
CommonLocs = cellfun(@(x) H.(x), CommonFields, 'uniformoutput', false)';
H = rmfield(H, CommonFields);
L = rmfield(L, CommonFields);

Cname = strrep(CommonFields, 'Loc', '')';
Hname = cellfun(@(x) ['h' strrep(x, 'Loc', '')], fieldnames(H), 'uniformoutput', false);
Lname = cellfun(@(x) ['l' strrep(x, 'Loc', '')], fieldnames(L), 'uniformoutput', false);

Bdata = [{'Chain'} Chain; Cname CommonLocs; Hname struct2cell(H); Lname struct2cell(L)]';
Map = struct(Bdata{:});

Map.ParNum = find(strcmpi(VDJheader, 'ParNum'));