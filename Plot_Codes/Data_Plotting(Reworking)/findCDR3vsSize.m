%
function PlotData = findCDR3vsSize()

[VDJdata,VDJheader,FileName,FilePath] = openSeqData;

H = getHeaderVar(VDJheader);

%Extract the unique groups
GrpNum = cell2mat(VDJdata(:,H.GrpNumLoc));
UnqGrpNum = unique(GrpNum);

Data = zeros(length(UnqGrpNum),3); %[GrpNum GrpSize CDR3Length]
for y = 1:length(UnqGrpNum)
    IdxLoc = find(UnqGrpNum(y) == GrpNum);
    Data(y,[1 3])= cell2mat(VDJdata(IdxLoc(1),[H.GrpNumLoc H.CDR3Loc(2)]));
    Data(y,2) = length(IdxLoc);
end

%Filter nonprod CDR3sw
DelThese = Data(:,3) == 0;
Data(DelThese,:) = [];

%Convert number of nts to number of aa
Data(:,3) = Data(:,3)/3;


X = 1:1:150;
Y = 3:35;
%Y = min(Data(:,3)):max(Data(:,3));
PlotData = zeros(length(X),length(Y));
Z = zeros(length(X),length(Y));
for x = 1:length(X)
    %Find how many clonotypes there are
    ValidClones = Data(:,2) >= X(x);
    Z(x,:) = sum(ValidClones);
    for y = 1:length(Y)
        PlotData(x,y) = sum(ValidClones & Data(:,3) == Y(y));
    end
end

NormPlotData = PlotData./Z;
imagesc(NormPlotData)

DotLoc = find(FileName == '.');
SaveName = [FileName(1:DotLoc(end)-1) '_CDR3vSize.png'];
saveas(gcf,[FilePath,SaveName]);
