%List of Common codes used over and over again in my other code
[VDJdata,VDJheader,FileName,FilePath] = openSeqData();

%Iteratively find groups
GrpNum = cell2mat(VDJdata(:,H.GrpNumLoc));
UnqGrpNum = unique(GrpNum);
for y = 1:length(UnqGrpNum)
    IdxLoc = find(UnqGrpNum(y) == GrpNum);
    
end

%Open multiple files
[FileNames,FilePath] = uigetfile('*.xlsx;*.csv','Open files','multiselect','on');
if ischar(FileNames)
    FileNames = {FileNames};
end
for f = 1:length(FileNames)
    [VDJdata,VDJheader,FileName,FilePath] = openSeqData([FilePath FileNames{f}]);
end

%Save the files
DotLoc = find(FileName == '.');
DotLoc = DotLoc(end);
SaveName = FileName(1:DotLoc-1);
%Before saving to xlsx, convert columns with matrix values into char
for q = 1:size(VDJdata,1)
    for w = 1:3
        VDJdata{q,H.FamNumLoc(w)} = mat2str(VDJdata{q,H.FamNumLoc(w)});
    end
end
if ispc
    xlswrite([FilePath SaveName 'D.xlsx'],[VDJheader; VDJdata]);
else
    writeDlmFile([VDJheader;VDJdata],[FilePath SaveName 'D.csv'],'\t');
end

%Input checking
if isempty(varargin)
    [VDJdata,VDJheader,FileName,FilePath] = openSeqData;
elseif length(varargin) == 2
    if ischar(varargin{1})
        FileName = varargin{1};
        FilePath = varargin{2};
        [VDJdata,VDJheader,FileName,FilePath] = openSeqData([FilePath FileName]);
    else
        VDJdata = varargin{1};
        VDJheader = varargin{2};
    end
end


%Extract the VDJ database
if length(varargin) >= 3
    [Vmap, Dmap, Jmap] = deal(varargin{1:3});
else
    [Vmap, Dmap, Jmap] = getCurrentDatabase;
end

%Parse the input
P = inputParser;
addRequired(P,'AncMap',@isnumeric);
addOptional(P,'TreeName','',@ischar);
addOptional(P,'CDR3seq','',@iscell);
addOptional(P,'DotClr',[],@isnumeric);
addParameter(P,'SortMode','sort',@(x) any(validatestring(x,{'sort','none'})));
addParameter(P,'FigWidth',3.3,@isnumeric);
parse(P,AncMap,varargin{:});
AncMap = P.Results.AncMap;
        
