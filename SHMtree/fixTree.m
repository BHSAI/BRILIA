%fixTree will look through VDJdata tree structures, and flip the root seq
%with first seq IF it doesn't make sense. One way to figure this out is to
%see if the 1st seq has MORE shm than another seq. Do this at the way end,
%after D and N region fixes.

function VDJdata = fixTree(VDJdata,NewHeader)
getHeaderVar;

GrpNum = cell2mat(VDJdata(:,GrpNumLoc));
UnqGrpNum = unique(GrpNum);
for y = 1:length(UnqGrpNum)
    IdxLoc = find(UnqGrpNum(y) == GrpNum);
       
    if length(IdxLoc) == 1; %Nothing to fix
        continue
    end

    Tdata = VDJdata(IdxLoc,:);
    SHMcount = sum(cell2mat(Tdata(:,SHMLoc)),2);
    MinLoc = find(SHMcount == min(SHMcount));
    TempCt = cell2mat(Tdata(MinLoc,TemplateLoc));
    CompareMat = [MinLoc TempCt];
    CompareMat = sortrows(CompareMat,-2);
    MinLoc = CompareMat(1,1);
    
    %if MinLoc == 1; continue; end %You want to avoid fixing tree if possible. The below linking is slighly worse than min-dist clustering method.
    %Turns out, redoing all tree improves SHM correlation. So just redo all
    %trees.
    
    %Change the root on AncMap
    %AncMapCell = calcAncMapCell(Tdata,NewHeader);
    AncMapCell = calcAncMapCell(Tdata,NewHeader);
    AncMap = AncMapCell{1};
    AncMap(1,2) = size(AncMap,1); %It can be anything, since we'll recalc this.
    AncMap(MinLoc,2) = 0;
       
    %Update the Root RefSeq.
    Tdata{MinLoc,RefSeqLoc} = Tdata{1,RefSeqLoc};
    AncMap = sortrows(AncMap,[2 1]);
    Tdata = Tdata(AncMap(:,1),:);
    
    %Update the refseqs after making the tree
    AncMap = buildTreeLink_B(Tdata,NewHeader);
    Tdata = mapData2AncMap(Tdata,NewHeader,AncMap,'SkipRef');
    
    VDJdata(IdxLoc,:) = Tdata;
end




