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
for j = 1:size(AncMap,1)
    ParLoc = find(PairDist(:,j) == min(PairDist(:,j)));
    AncMap(j,1:3) = [j ParLoc(1) min(PairDist(:,j))];
end
AncCycle = findTreeCycle(AncMap);


ChangedClust = 1; %Counter to detect changes to clusteringj
KillTimer2 = 1;
while max(AncCycle) > 0
    if KillTimer2 > 10000
        disp('KillTimer2');
        break
    else
        KillTimer2 = KillTimer2+1;
    end
    
    %Find the root seq per cluster, breaking up cyclic dependencies.
    AncMap(AncMap(:,3) > CutoffDist,2) = 0; %Deactivating large distances
    ClustMap = findTreeClust(AncMap); 
    AncMap(:,5) = ClustMap(:,2);
    AncMap = sortrows(AncMap,[5 2]);

    for y = 1:max(AncMap(:,5))
        IdxLoc = find(AncMap(:,5) == y);
        AncMapT = AncMap(IdxLoc,:);
        CycleLoc = findTreeCycle(AncMapT); %Find the cyclic dependency and also minimum distance

        if max(CycleLoc) > 0 %This means you need to find root
            %Root is the sequence with the smallest distance to other seqs
            PairDistT = PairDist(AncMapT(:,1),:);
            PairDistT = PairDistT(:,AncMapT(:,1));
            PairDistT(eye(size(PairDistT))>0) = 0;
            TotDist = sum(PairDistT,2);
            RootLoc = AncMapT(TotDist == min(TotDist(CycleLoc)),1);
            
            %If there are multiple roots, choose shortest SHM dist to RefSeq, and then largest template count.
            if length(RootLoc) > 1           
                VDJscore = zeros(length(RootLoc),5);
                for g = 1:length(RootLoc)
                    %Extract data needed to calculate scores
                    VMDNJ = cell2mat(Tdata(RootLoc(g),LengthLoc));
                    CurSeq = Tdata{RootLoc(g),SeqLoc};
                    RefSeq = Tdata{RootLoc(g),RefSeqLoc};
                    SeqDiff = CurSeq == RefSeq;
                    TempCt = AncMap(RootLoc(g)==AncMap(:,1),4);
                    
                    %Compute the scores
                    Vscore = calcAlignScore(SeqDiff(1:VMDNJ(1)),inf);
                    Dscore = calcAlignScore(SeqDiff(sum(VMDNJ(1:2))+1:sum(VMDNJ(1:3))),inf);
                    Jscore = calcAlignScore(SeqDiff(sum(VMDNJ(1:4))+1:sum(VMDNJ)),inf);

                    %Fill in the score matrix
                    VDJscore(g,:) = [Vscore Jscore Dscore TempCt  RootLoc(g)];
                end
                VDJscore = sortrows(VDJscore,[-1 -2 -3 -4]); %Sort by highest alignment scores, then largest template count.
                RootLoc = VDJscore(1,end); %Get the first one for now   
            end
            AncMapT(AncMapT(:,1) == RootLoc,2) = 0; %Set parent to 0
        end
        
        %Update the main table and counter
        AncMap(IdxLoc,:) = AncMapT;
    end
    AncMap = sortrows(AncMap,[5 2]);

    %----------------------------------------------------------------------
    %All subclusters should have roots now. Need to cluster-cluster link
    RootLocs = AncMap(AncMap(:,2) == 0,1);
    PairDistC = zeros(length(RootLocs)); %Calculate cluster-cluster distances
    NodeLocC = zeros(length(RootLocs)); %Keeps track of which node is the closest to the root
    for r = 1:size(PairDistC,1)
        for c = 1:size(PairDistC,2)
            if r == c; continue; end
            PosParent = AncMap(AncMap(:,5) == r,1); %Possible parent of this cluster c
            PairDistT = PairDist(PosParent,RootLocs(c)); %Extract distance of all potential paren to this cluster c
            NodeLocs = find(PairDistT == min(PairDistT(:))); %Find the closest possible parent
            NodeLocC(r,c) = PosParent(NodeLocs(1)); %Save closest possible parent
            PairDistC(r,c) = PairDistT(NodeLocs(1)); %Save distance to closest possible parent
        end
    end
    PairDistC(eye(size(PairDistC))>0) = inf; %Prevent self match

    %Cluster the clusters
    AncMapC = zeros(size(PairDistC,1),5); %[ChildNum ParNum ChildtoParNodeDist ParNodeNum ClustNum]
    for j = 1:size(AncMapC,1)
        ParLocC = find(PairDistC(:,j) == min(PairDistC(:,j)));
        AncMapC(j,:) = [j ParLocC(1) min(PairDistC(:,j)) NodeLocC(ParLocC(1),j) 0];
    end
    AncMapC(AncMapC(:,3) > CutoffDist,2) = 0; %Cutoff distance
    ClustMapC = findTreeClust(AncMapC);
    AncMapC(:,end) = ClustMapC(:,2);
    
    %For all linkable clusters, consolidate group number, fill in the
    %parent seq. Let the previous for loop above take care of root search.
    for w = 1:max(AncMapC(:,end))
        ClustTLoc = AncMapC(:,end)==w;
        CurRoot = RootLocs(ClustTLoc);
        ParNode = AncMapC(ClustTLoc,4);
        ParDist = AncMapC(ClustTLoc,3);
        if length(CurRoot) >  1
            for q = 1:length(CurRoot)
                AncMap(AncMap(:,1) == CurRoot(q),2:3) = [ParNode(q) ParDist(q)]; %Change this cluster's root 0 to node of other cluster. Reset child2par distance too.
            end
        end
    end
    AncCycle = findTreeCycle(AncMap);
    
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
    
%     AncMapCell = calcAncMapCell(TdataNew,NewHeader);
%     for k = 1:size(AncMapCell,1)
%         NewAncMap = AncMapCell{k};
%         TreeMapTemp= findTreeClust(NewAncMap);
%         if max(TreeMapTemp(:,2)) > 1; 
%             error('You have multiple clusters')
%         elseif sum(NewAncMap(:,1) == NewAncMap(:,2)) > 0
%             error('parent childe same')
%         elseif sum(NewAncMap(:,2) == 0 ) > 1
%             error('multiple Refs')
%         elseif sum(NewAncMap(:,3) == 0 ) > 1
%             error('multiple 0 dist')
%         end
%     end
%     
    varargout{1} = TdataNew;
end