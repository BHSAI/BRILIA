%renumberGrpNum will renumber the group # to go from 1 to N.
%
%  VDJdata = renumberGrpNum()
%
%  VDJdata = renumberGrpNum(FullFileName)
%
%  VDJdata = renumberGrpNum(VDJdata,VDJheader);
%

function VDJdata = renumberGrpNum(varargin)
%Prepare the inputs
if isempty(varargin)
    [VDJdata,VDJheader,FileName,FilePath] = openSeqData;
else
    if ischar(varargin{1})
        FullFileName = varargin{1};
        [VDJdata,VDJheader,FileName,FilePath] = openSeqData(FullFileName);
    elseif iscell(varargin{1})
        VDJdata = varargin{1};
        VDJheader = varargin{2};
        FileName = '';
        FilePath = '';
    end
end
H = getHeaderVar(VDJheader);

%Begin renumbering the group numbers
GrpCt = 1;
CurGrp = VDJdata{1,H.GrpNumLoc};
for j = 1:size(VDJdata,1)
    if VDJdata{j,H.GrpNumLoc} ~= CurGrp
        GrpCt = GrpCt + 1;
        CurGrp = VDJdata{j,H.GrpNumLoc};
    end
    VDJdata{j,H.GrpNumLoc} = GrpCt;
end

%Save the file
if ~isempty(FileName)
    DotLoc = find(FileName == '.');
    if isempty(DotLoc)
        DotLoc = length(FileName);
    end
    SaveDelimiter = ';';
    SaveFullName = sprintf('%s%s.Renum.%s',FilePath,FileName(1:DotLoc(end)-1),'csv');
    saveSeqData(SaveFullName,VDJdata,VDJheader,'Delimiter',SaveDelimiter);
end
