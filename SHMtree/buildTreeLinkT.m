%buildAncMap will take a cluster of sequences, and assemble a tree such
%that the maximum distance of X is preserved according to the
%calcSHMHAMdist function.

%1) Cluster based on distance first, cyclic dependencies allowed

%2) For each cluster, identify the "root sequence". This is the sequences
%that gives you the smallest total distance to all other sequences.

%3) For each root of each cluster, attempt to link to another cluster's
%node or leaf.

function [AncMapCell,varargout] = buildTreeLink(Tdata,NewHeader,varargin)
%Remember, pairwise distance from SHMHAM is doubled the hamming distance
%because original SHMHAM dist uses 0.5 fraction, which is bad for memory
%usage of matrices that could have used int16 format. 
if isempty(varargin)
    CutoffDist = 4*2; 
else
    CutoffDist = varargin{1}*2;
end
getHeaderVar;

%==========================================================================
%Identify starting point clusters. Will do cluster-cluster linking later.
[PairDist,~] = calcPairDist(Tdata(:,SeqLoc),'shmham'); %Parent is each rows, Child is each column. Note that SHMHAM distance is doubled by default, to ensure integer values.
AncMap = [[1:size(PairDist,1)]' zeros(size(PairDist,1),4)]; %[ChildNum ParNum SHMHAMdist Template ClusterNum]
if TemplateLoc > 0
    AncMap(:,4) = cell2mat(Tdata(:,TemplateLoc));
    AncMap(isnan(AncMap(:,4)),4) = 1;
    AncMap(AncMap(:,4)==0,4) = 1;
else
    AncMap(:,4) = 1;
end

%Link AncMap by nearest distance neighbor (cyclic dependencies inevitable)
for j = 1:size(AncMap,1)
    ParLoc = find(PairDist(:,j) == min(PairDist(:,j)));
    AncMap(j,1:3) = [j ParLoc(1) min(PairDist(:,j))];
end

%Find the tree clusters 
AncMap(AncMap(:,3) > CutoffDist,2) = 0; %Deactivating large distances
ClustMap = findTreeClust(AncMap); 
AncMap(:,5) = ClustMap(:,2);

%For each cluster, identify the "root" as one with V>J>D scores and tempct
y = 1;
while y <= max(AncMap(:,5))
    AncIdx = find(AncMap(:,5) == y);
    AncMapT = AncMap(AncIdx,:);
    DataIdx = AncMap(AncIdx,1);
    VDJTscore = zeros(length(DataIdx),5);
    for g = 1:length(AncIdx)
        %Extract data needed to calculate scores
        VMDNJ = cell2mat(Tdata(DataIdx(g),LengthLoc));
        CurSeq = Tdata{DataIdx(g),SeqLoc};
        RefSeq = Tdata{DataIdx(g),RefSeqLoc};
        SeqDiff = CurSeq == RefSeq;
        TempCt = AncMapT(g,4);
                    
        %Compute the scores
        Vscore = calcAlignScore(SeqDiff(1:VMDNJ(1)),inf);
        Dscore = calcAlignScore(SeqDiff(sum(VMDNJ(1:2))+1:sum(VMDNJ(1:3))),inf);
        Jscore = calcAlignScore(SeqDiff(sum(VMDNJ(1:4))+1:sum(VMDNJ)),inf);

        %Fill in the score matrix
        VDJTscore(g,:) = [Vscore Jscore Dscore TempCt g];
    end
    VDJTscore = sortrows(VDJTscore,[-1 -2 -3 -4]); %Sort by highest alignment scores, then largest template count.
    RootAncIdx = VDJTscore(1,end); %Get the first one for now       
    AncMapT(RootAncIdx,2) = 0; %Set this to root;
    AncMapT([1 RootAncIdx],:) = AncMapT([RootAncIdx 1],:); %Just make sure first one is root.
    
    %If setting the root creates more clusters, renumber to new clusters
    AncClustT = findTreeClust(AncMapT);
    AncClustT(AncClustT(:,end)>1,end) = AncClustT(AncClustT(:,end)>1,end) + max(AncMap(:,5)) - 1;
    AncClustT(AncClustT(:,end)==1,end) = y;
    AncMapT(:,end) = AncClustT(:,end);
    
    %Fill back the AncMap info
    AncMap(AncIdx,:) = AncMapT;
    y = y+1;
end
AncMap = sortrows(AncMap,[5 2]);

%----------------------------------------------------------------------
%All subclusters should have roots now. Need to cluster-cluster link
NewClust = 1;
KillTimer = 1;
while NewClust
    KillTimer = KillTimer+1;
    if KillTimer > 100
        disp('KillTimer breaking while loop')
        pause
        break
    end

    RootLocs = AncMap(AncMap(:,2) == 0,1);
    
    %If there is only 1 root, then just stop here
    if length(RootLocs) == 1; break; end
    
    PairDistC = zeros(length(RootLocs)); %Calculate cluster-cluster distances
    NodeLocC = zeros(length(RootLocs)); %Keeps track of which node is the closest to the root
    for r = 1:size(PairDistC,1)
        for c = 1:size(PairDistC,2)
            if r == c; continue; end
            PosParent = AncMap(AncMap(:,end) == r,1); %Possible parent of this cluster c
            PairDistT = PairDist(PosParent,RootLocs(c)); %Extract distance of all potential paren to this cluster c
            NodeLocs = find(PairDistT == min(PairDistT(:))); %Find the closest possible parent
            NodeLocC(r,c) = PosParent(NodeLocs(1)); %Save closest possible parent
            PairDistC(r,c) = PairDistT(NodeLocs(1)); %Save distance to closest possible parent
        end
    end
    PairDistC(eye(size(PairDistC))>0) = inf; %Prevent self match

    %Cluster the clusters
    AncMapC = zeros(size(PairDistC,1),6); %[ChildNum ParNum ChildtoParNodeDist ParNodeNum GermRefRootDist ClustNum]
    for j = 1:size(AncMapC,1)
        ParLocC = find(PairDistC(:,j) == min(PairDistC(:,j)));

        %Find each cluster's root SHMdist to germline
        RootSeqLoc = find(AncMap(:,5) == j);
        RootSeqIdx = AncMap(RootSeqLoc(1),1);
        GermSeq = Tdata{RootSeqIdx,RefSeqLoc};
        RootSeq = Tdata{RootSeqIdx,SeqLoc};
        GermRootDist = calcSHMHAMdist(GermSeq,RootSeq);

        AncMapC(j,:) = [j ParLocC(1) min(PairDistC(:,j)) NodeLocC(ParLocC(1),j) GermRootDist 0];
    end
    AncMapC(AncMapC(:,3) > CutoffDist,2) = 0; %Cutoff distance
    ClustMapC = findTreeClust(AncMapC); 
    AncMapC(:,end) = ClustMapC(:,2);
    AncMapC = sortrows(AncMapC,[6 2]);

    %For every cyclic dependency, break the furthest distance one. 
    for t = 1:max(AncMapC(:,end))
        ClustIdx = find(AncMapC(:,end) == t);
        AncMapCT = AncMapC(ClustIdx,:);
        CycleLoc = find(findTreeCycle(AncMapCT));
        if ~isempty(CycleLoc)
            RootClust = find(AncMapCT(CycleLoc,2) == max(AncMapCT(CycleLoc,5))); %Greater distance loses.
            AncMapCT(CycleLoc(RootClust(1)),2) = 0;
            AncMapC(ClustIdx,:) = AncMapCT;
        end
    end

    %Link the AncMap together based on AncMapC
    NewClust = 0;
    for k = 1:size(AncMapC,1)
        if AncMapC(k,2) > 0
            NewClust = 1;
            AncIdx = find(AncMap(:,end) == k);
            AncMap(AncIdx(1),2:3) = AncMapC(k,[4 3]);
        end
    end
end

%Perform final sorting and cutting
AncMap(AncMap(:,3) > CutoffDist,2) = 0; %Deactivating large distances
ClustMap = findTreeClust(AncMap); 
AncMap(:,5) = ClustMap(:,2);
AncMap = sortrows(AncMap,[5 2]);

%Basic formatting
UnqClust = unique(AncMap(:,end));
AncMapCell = cell(length(UnqClust),2);
for j= 1:size(AncMapCell,1)
    AncMapT = AncMap(AncMap(:,end)==j,:);
    AncMapT(:,3) = AncMapT(:,3)/2; %Remember, SHMHAM was doubled for algorithm purposes, so need to divide back by 2;
    AncMapCell{j,1} = AncMapT;
    AncMapCell{j,2} = Tdata(AncMapT(:,1),CDR3Loc(1)); %Save the CDR3 info
end

%Restructure Tdata according to ordering
if nargout == 2   
    GrpNum = 1;
    S1 = 1; %Start index of TdataNew
    TdataNew = cell(size(Tdata,1),length(NewHeader));
    for k = 1:size(AncMapCell,1)
        
        %Extract the data for this cluster
        AncMapT = AncMapCell{k,1};
        TdataT = Tdata(AncMapT(:,1),:);
        
        %Renumber each AncMapCell so that AncMapT(:,1) is 1,2,3...
        for h = 1:size(AncMapT,1)
            ParLoc = find(AncMapT(h,2) == AncMapT(:,1));
            if ~isempty(ParLoc)
                AncMapT(h,2) = ParLoc;
            end
        end
        AncMapT(:,1) = 1:size(AncMapT,1);
        AncMapCell{k} = AncMapT;
        
        %If there is the TreeCountLoc, fill in the data here
        if ChildCountLoc > 0
            for w = 1:size(TdataT,1)
                TdataT{w,ChildCountLoc} = length(findChild(AncMapT,AncMapT(w,1)));
            end
        end
            
        %Switch off the RefSeq with the ParentSeq
        for h = 1:size(AncMapT,1)
            CurLoc = AncMapT(h,2);
            if CurLoc > 0
                ParLoc = AncMapT(:,1) == CurLoc;
                TdataT(h,RefSeqLoc) = TdataT(ParLoc,SeqLoc);
            end
        end
        
        %Add the TdataNew entries, and update S1 and GrpNum counters
        S2 = S1 + size(AncMapCell{k,1},1)-1; %Final index of TdataNew
        TdataNew(S1:S2,1:size(TdataT,2)) = TdataT;
        TdataNew(S1:S2,GrpNumLoc) = repmat({GrpNum},S2-S1+1,1);
        S1 = S2 + 1;
        GrpNum = GrpNum + 1;
    end
    
    AncMapCell = calcAncMapCell(TdataNew,NewHeader);
    for k = 1:size(AncMapCell,1)
        NewAncMap = AncMapCell{k};
        TreeMapTemp= findTreeClust(NewAncMap);
        if max(TreeMapTemp(:,2)) > 1; 
            error('You have multiple clusters')
        elseif sum(NewAncMap(:,1) == NewAncMap(:,2)) > 0
            error('parent childe same')
        elseif sum(NewAncMap(:,2) == 0 ) > 1
            error('multiple Refs')
        elseif sum(NewAncMap(:,3) == 0 ) > 1
            error('multiple 0 dist')
        end
    end
    
    varargout{1} = TdataNew;
end


