%renumberGrpNum will renumber the group # to go from 1 to N.
%
%  VDJdata = renumberGrpNum()
%
%  VDJdata = renumberGrpNum(FullFileName)
%
%  VDJdata = renumberGrpNum(VDJdata,NewHeader);
%

function VDJdata = renumberGrpNum(varargin)
%Prepare the inputs
if isempty(varargin)
    [VDJdata,NewHeader,FileName,FilePath] = openSeqData;
else
    if ischar(varargin{1})
        FullFileName = varargin{1};
        [VDJdata,NewHeader,FileName,FilePath] = openSeqData(FullFileName);
    elseif iscell(varargin{1})
        VDJdata = varargin{1};
        NewHeader = varargin{2};
        FileName = '';
        FilePath = '';
    end
end
getHeaderVar;

%Begin renumbering the group numbers
GrpCt = 1;
CurGrp = VDJdata{1,GrpNumLoc};
for j = 1:size(VDJdata,1)
    if VDJdata{j,GrpNumLoc} ~= CurGrp
        GrpCt = GrpCt + 1;
        CurGrp = VDJdata{j,GrpNumLoc};
    end
    VDJdata{j,GrpNumLoc} = GrpCt;
end

%Save the file
if ~isempty(FileName)
    DotLoc = find(FileName == '.');
    if isempty(DotLoc)
        DotLoc = length(FileName);
    end
    SaveDelimiter = ';';
    SaveFullName = sprintf('%s%s.Renum.%s',FilePath,FileName(1:DotLoc(end)-1),'csv');
    saveSeqData(SaveFullName,VDJdata,NewHeader,'Delimiter',SaveDelimiter);
end
