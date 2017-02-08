%getTreeStats will return the following tree-related info: [GroupNum
%WtNodeLevel AvgChildperNode]

function TreeStats = findTreeStats(varargin)
if isempty(varargin)
    [VDJdata,VDJheader,FileName,FilePath] = openSeqData;
else
    VDJdata = varargin{1};
    VDJheader = varargin{2};
end
H = getHeaderVar(VDJheader);

GrpNum = cell2mat(VDJdata(:,H.GrpNumLoc));
GrpNumUnq = unique(GrpNum);
TreeStats = zeros(length(GrpNumUnq),7); %ClustNum | WtNodeLevel | AvgChildperNode | MaxMutFullBCR| NumSeq | MaxHamDistTraveled | NodeCount
for y = 1:length(GrpNumUnq)
    IdxLoc = find(GrpNumUnq(y) == GrpNum);    
    if length(IdxLoc) == 1; %No point in roots a single sequences
        continue
    end
    Tdata = VDJdata(IdxLoc,:);
    DistMode = 'shmham';
    [AncMapCell,~] = buildTreeLink3(Tdata,VDJheader,DistMode);

    AncMap = AncMapCell{1,1};
    if AncMap(1,2) ~= 0
        disp('Setting 1st seq to be linked to root, Parent 0');
        AncMap(1,2) = 0;
    end
    
    NodeData = zeros(size(AncMap,1),2); %[NodeLevel ChildCt]
    for j = 1:size(AncMap,1)
        NodeData(j,1) = findSeqNodeLevel(AncMap,AncMap(j,1));
        NodeData(j,2) = length(findChild(AncMap,AncMap(j,1)));
    end
    
    RefSeq = Tdata{1,H.RefSeqLoc};
    SamSeq = char(Tdata(:,H.SeqLoc));
    
    %Compute maximum hamming distance to RefSeq
    CmprRefSamSeq = SamSeq ~= repmat(RefSeq,size(SamSeq,1),1);
    HamDists = sum(CmprRefSamSeq,2);
    MaxMutations = max(HamDists(:))/size(RefSeq,2);
    
    %COmpute the maximum hamming distance traveled from germline
    TreeCoord = calcTreeCoord(AncMap);
    
    
    TreeStats(y,1) = GrpNumUnq(y);
    TreeStats(y,2) = max(NodeData(:,1)); %Max NodeLevel
    TreeStats(y,3) = max(NodeData(:,2)); %Max Children per node
    TreeStats(y,4) = MaxMutations; %Maximum extent of mutation
    TreeStats(y,5) = size(Tdata,1); %Number of unique sequences
    TreeStats(y,6) = max(TreeCoord(:,1)); %Number of unique sequences
    TreeStats(y,7) = sum(TreeCoord(:,3) == 0); %Number of unique sequences
end

DelThese = TreeStats(:,1) == 0;
TreeStats(DelThese,:) = [];

if exist('FileName','var')

    DotLoc = find(FileName == '.');
    FilePre = FileName(1:DotLoc-1);
    save([FilePath FilePre 'Tree.mat'],'TreeStats')
end
