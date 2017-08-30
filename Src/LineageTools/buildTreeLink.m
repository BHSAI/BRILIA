%buildTreeLink will take a crude cluster of VDJdata (Tdata) and then use
%lineage trees to determine the fine clusters based on the SHM distance.
%This will return a rearranged Tdata cell matrix based on the new
%clustering, where the 1st member is always the sequence closest to the
%germline sequence.
%
%Procedure are as follows:
%1) Cluster based on distance first, cyclic dependencies allowed
%
%2) For each cluster, identify the "root sequence". This is the sequences
%that gives you the smallest total SHM distance to all other sequences.
%
%3) For each root of each cluster, attempt to link to another cluster's
%node or leaf. Repeat until clusters do not change.
%
%  [AncMapCell,Tdata] = buildTreeLink(Tdata,VDJheader,DevPerc)
%
%  AncMapCell = buildTreeLink(Tdata,VDJheader,DevPerc)
%
%  INPUT
%    Tdata: main BRILIA data cell, but selecting for only sequences that
%      belong in the same crude cluster (see clusterGene.m) and have the
%      same sequence lengths (padtrimSeqGroup.m);
%    VDJheader: main BRILIA header cell
%    DevPerc: Ranges 0 to 100, indicating % of sequence length to use as
%      the cutoff SHM distance. 
%    S: special header structure containing locations of important data
%      columns. This is used to speed up the repetitive header search. Can
%      get this using S = buildTreeLink('getHeaderVar',VDJheader);
%        S.TemplateLoc
%        S.TreeChildLoc 
%        S.HSeqLoc
%        S.LSeqLoc 
%        S.HRefSeqLoc 
%        S.LRefSeqLoc 
%
%  OUTPUT
%    AncMapCell: a Mx1 cell, where each cell contain an ancestry map matrix
%      (AncMap) that containing the following information:
%        [ChildSeqNum  ParentSeqNum  HAMdist  SHMdist]
%    Tdata: rearranged Tdata such that the clusters are grouped together,
%      with the 1st sequence of each cluster being closest to the germline
%      seq.


function varargout = buildTreeLink(varargin)
%See if user wants to parse the header (special function)
JustGettingS = 0;
S = [];
if nargin >= 2 && ischar(varargin{1}) && strcmpi(varargin{1},'getheadervar')
    VDJheader = varargin{2};
    JustGettingS = 1;
elseif nargin >= 3
    Tdata = varargin{1};
    VDJheader = varargin{2};
    DevPerc = varargin{3};
    if nargin >= 4 && isstruct(varargin{4})
        S = varargin{4};
    end
else
    error('%s: Wrong number of inputs',mfilename);    
end

%Extract the minimally needed header information due to high repetetition
if isempty(S)
    H = getHeavyHeaderVar(VDJheader);
    L = getLightHeaderVar(VDJheader);
    S.TemplateLoc = H.TemplateLoc;
    S.ChildCountLoc = H.ChildCountLoc;
    S.GrpNumLoc = H.GrpNumLoc;
    S.HCDR3Loc = H.CDR3Loc;
    S.LCDR3Loc = L.CDR3Loc;
    S.HSeqLoc = H.SeqLoc;
    S.LSeqLoc = L.SeqLoc;
    S.HRefSeqLoc = H.RefSeqLoc;
    S.LRefSeqLoc = L.RefSeqLoc;
    S.HLengthLoc = H.LengthLoc;
    S.LLengthLoc = L.LengthLoc;
end

%Return S if specified
if JustGettingS
    varargout{1} = S;
    return
end

%Determine the chain
if S.HSeqLoc > 0 && S.LSeqLoc > 0
    Chain = 'HL';
elseif S.HSeqLoc > 0
    Chain = 'H';
else
    Chain = 'L';
end

%Determine the maximum Seq Dist, which is used to set the class of PairDist
MaxDist = 0;
MaxSeqLen = 0;
for c = 1:length(Chain)
    MaxDist = MaxDist + length(Tdata{1,S.([Chain(c) 'SeqLoc'])}) ^ 2+length(Tdata{1,S.([Chain(c) 'SeqLoc'])});
    MaxSeqLen = MaxSeqLen + length(Tdata{1,S.([Chain(c) 'SeqLoc'])});
end
Class = 'uint16';
if MaxDist > 2^32-1
    Class = 'double';
end

%Determine the distances between sequences
PairDist = zeros(size(Tdata,1),Class);
for c = 1:length(Chain) %Add H and L chain distances
    PairDist = PairDist + calcPairDist(Tdata(:,S.([Chain(c) 'SeqLoc'])),'shmham',Class); %Seq-Seq distances
end

%Compute the cutoff distance
CutoffDist = ceil(MaxSeqLen*DevPerc/100) * 2; %Remember that the PairDist is DOUBLED!!!

%Calculate the cutoff SHM distance, accounting for the fact that we must
%double the value because pairwise SHM distance is doubled for memory
%limitations. The original SHM dist in the paper uses 0.5 fractions, which
%would have required the use of a double matrix. Doubling the value allows
%for the use of uint matrices instead.

%==========================================================================
% %Ensure all sequence lengths are the same, aligned based on the CDR3endLoc.
% Tdata = trimSeq(Tdata,VDJheader);

%Determine initial parent-child pairing
if S.TemplateLoc > 0
    TempCt = cell2mat(Tdata(:,S.TemplateLoc));
else
    TempCt = ones(size(Tdata,1),1);    
end
ClustNum = zeros(size(Tdata,1),1);
AncMap = [calcAncMap(PairDist) TempCt ClustNum];  %Nearest neighbor method to links 
AncCycle = findTreeCycle(AncMap);
KillTimer2 = 1;
while max(AncCycle) > 0
    %fprintf('%s: AncCycle left %d\n',mfilename,sum(AncCycle));
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
        IdxLoc = AncMap(:,5) == y;
        AncMapT = AncMap(IdxLoc,:);
        CycleLoc = findTreeCycle(AncMapT); %Find the cyclic dependency and also minimum distance

        if max(CycleLoc) > 0 %This means you need to find root
            %Root is the sequence with the smallest distance to other seqs
            PairDistT = PairDist(AncMapT(:,1),AncMapT(:,1)); 
            TotDist = sum(PairDistT,2);
            RootLoc = AncMapT(TotDist == min(TotDist(CycleLoc)),1);
            
            %If there are multiple roots, choose shortest SHM dist to RefSeq, and then largest template count.
            if length(RootLoc) > 1
                if ~isempty(strfind(Chain,'H'))
                    VDJscore = zeros(length(RootLoc),5);
                    for g = 1:length(RootLoc)
                        %Extract data needed to calculate scores
                        VMDNJ = cell2mat(Tdata(RootLoc(g),S.HLengthLoc));
                        CurSeq = Tdata{RootLoc(g),S.HSeqLoc};
                        RefSeq = Tdata{RootLoc(g),S.HRefSeqLoc};
                        
                        %Make sure all info is there
                        if isempty(VMDNJ); continue; end
                        if isempty(CurSeq); continue; end
                        if isempty(RefSeq); continue; end
                        if length(CurSeq) ~= length(RefSeq); continue; end
                        
                        SeqDiff = CurSeq == RefSeq;

                        %Compute the scores
                        Vscore = calcAlignScore(SeqDiff(1:VMDNJ(1)),inf);
                        Dscore = calcAlignScore(SeqDiff(sum(VMDNJ(1:2))+1:sum(VMDNJ(1:3))),inf);
                        Jscore = calcAlignScore(SeqDiff(sum(VMDNJ(1:4))+1:sum(VMDNJ)),inf);

                        %Fill in the score matrix
                        TempCt = AncMap(RootLoc(g)==AncMap(:,1),4);
                        VDJscore(g,:) = [Vscore Dscore Jscore TempCt RootLoc(g)];
                    end
                end
                if ~isempty(strfind(Chain,'L'))
                    VJscore = zeros(length(RootLoc),4);
                    for g = 1:length(RootLoc)
                        %Extract data needed to calculate scores
                        VNJ = cell2mat(Tdata(RootLoc(g),S.LLengthLoc));
                        CurSeq = Tdata{RootLoc(g),S.LSeqLoc};
                        RefSeq = Tdata{RootLoc(g),S.LRefSeqLoc};

                        %Make sure all info is there
                        if isempty(VNJ); continue; end
                        if isempty(CurSeq); continue; end
                        if isempty(RefSeq); continue; end
                        if length(CurSeq) ~= length(RefSeq); continue; end

                        SeqDiff = CurSeq == RefSeq;

                        %Compute the scores
                        Vscore = calcAlignScore(SeqDiff(1:VNJ(1)),Inf);
                        Jscore = calcAlignScore(SeqDiff(sum(VNJ(1:2))+1:sum(VNJ)),Inf);

                        %Fill in the score matrix
                        TempCt = AncMap(RootLoc(g)==AncMap(:,1),4);
                        VJscore(g,:) = [Vscore Jscore TempCt RootLoc(g)];
                    end                    
                end
                if strcmpi(Chain,'HL')
                    AllScore = [VDJscore(:,1:3) VJscore];
                    RankOrder = -[1 3 4 5 2 6]; %V>J>Vx>Jx>D>TempCt
                elseif strcmpi(Chain,'H')
                    AllScore = VDJscore;
                    RankOrder = -[1 3 2 4]; %V>J>D>TempCt
                else
                    AllScore = VJscore;
                    RankOrder = -[1 2 3]; %Vx>Jx>TempCt
                end
                
                AllScore = sortrows(AllScore,RankOrder); %Sort by highest alignment scores, then largest template count.
                RootLoc = AllScore(1,end); %Get the first one for now\
                if RootLoc == 0; RootLoc = 1; end %This shouldn't happen, but if it does, proceed with 1st one.
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
        PosParent = AncMap(AncMap(:,5) == r,1); %Possible parent of this cluster c
        for c = 1:size(PairDistC,2)
            if r == c; continue; end
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
AncMapCell = cell(length(UnqClust),3);
for j= 1:size(AncMapCell,1)
    AncMapT = AncMap(AncMap(:,end)==j,:);
    AncMapT(:,3) = AncMapT(:,3)/2; %Remember, SHMHAM was doubled for algorithm purposes, so need to divide back by 2;
    AncMapCell{j,1} = AncMapT;
    if ~isempty(strfind(Chain,'H'))
        AncMapCell{j,2} = Tdata(AncMapT(:,1),S.HCDR3Loc(1)); %Save the CDR3 info
    end
    if ~isempty(strfind(Chain,'L'))
        AncMapCell{j,3} = Tdata(AncMapT(:,1),S.LCDR3Loc(1)); %Save the CDR3 info
    end
end

%Restructure Tdata according to ordering
if nargout >= 1
    varargout{1} = AncMapCell;
    if nargout == 2   
        GrpNum = 1;
        S1 = 1; %Start index of TdataNew
        TdataNew = cell(size(Tdata,1),length(VDJheader));
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
            if S.ChildCountLoc > 0
                for w = 1:size(TdataT,1)
                    TdataT{w,S.ChildCountLoc} = length(findChild(AncMapT,AncMapT(w,1)));
                end
            end

            %Switch the RefSeq with the ParentSeq
            for h = 1:size(AncMapT,1)
                CurLoc = AncMapT(h,2);
                if CurLoc > 0
                    for c = 1:length(Chain)
                        ParLoc = AncMapT(:,1) == CurLoc;
                        TdataT(h,S.([Chain(c) 'RefSeqLoc'])) = TdataT(ParLoc,S.([Chain(c) 'SeqLoc']));
                    end
                end
            end

            %Add the TdataNew entries, and update S1 and GrpNum counters
            S2 = S1 + size(AncMapCell{k,1},1)-1; %Final index of TdataNew
            TdataNew(S1:S2,1:size(TdataT,2)) = TdataT;
            TdataNew(S1:S2,S.GrpNumLoc) = repmat({GrpNum},S2-S1+1,1);
            S1 = S2 + 1;
            GrpNum = GrpNum + 1;
        end
        varargout{2} = TdataNew;
    end
end
