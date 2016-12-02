%renumberGrpNum will take VDJdata, and renumber the groups to go from 1 to
%N.

function VDJdata = renumberGrpNum(varargin)
if isempty(varargin)
    [VDJdata,NewHeader,FileName,FilePath] = openSeqData;
else
    VDJdata = varargin{1};
    NewHeader = varargin{2};
    FileName = '';
    FilePath = '';
end
getHeaderVar;

GrpCt = 1;
CurGrp = VDJdata{1,GrpNumLoc};
for j = 1:size(VDJdata,1)
    if VDJdata{j,GrpNumLoc} ~= CurGrp
        GrpCt = GrpCt + 1;
        CurGrp = VDJdata{j,GrpNumLoc};
    end
    VDJdata{j,GrpNumLoc} = GrpCt;
end

if ~isempty(FileName)
    %Save the files
    DotLoc = find(FileName == '.');
    DotLoc = DotLoc(end);
    SaveName = FileName(1:DotLoc-1);
    %Before saving to xlsx, convert columns with matrix values into char
    for q = 1:size(VDJdata,1)
        for w = 1:3
            VDJdata{q,FamNumLoc(w)} = mat2str(VDJdata{q,FamNumLoc(w)});
        end
    end
    if ispc
        xlswrite([FilePath SaveName 'D.xlsx'],[NewHeader; VDJdata]);
    else
        writeDlmFile([NewHeader;VDJdata],[FilePath SaveName 'D.csv'],'\t');
    end
end