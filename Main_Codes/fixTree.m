%fixTree will look through VDJdata tree structures, and flip the root seq
%with first seq IF it doesn't make sense. One way to figure this out is to
%see if the 1st seq has MORE shm than another seq. Do this at the way end,
%after D and N regions are fixed.
%
%  VDJdata = fixTree(VDJdata,VDJheader)

function VDJdata = fixTree(VDJdata,VDJheader)
H = getHeaderVar(VDJheader);

GrpNum = cell2mat(VDJdata(:,H.GrpNumLoc));
UnqGrpNum = unique(GrpNum);
for y = 1:length(UnqGrpNum)
    try
        IdxLoc = find(UnqGrpNum(y) == GrpNum);

        if length(IdxLoc) == 1; %Nothing to fix
            continue
        end

        Tdata = VDJdata(IdxLoc,:);
        SHMcount = sum(cell2mat(Tdata(:,H.SHMLoc)),2);
        MinLoc = find(SHMcount == min(SHMcount));
        TempCt = cell2mat(Tdata(MinLoc,H.TemplateLoc));
        CompareMat = [MinLoc TempCt];
        CompareMat = sortrows(CompareMat,-2);
        MinLoc = CompareMat(1,1);

        %if MinLoc == 1; continue; end %You want to avoid fixing tree if possible. The below linking is slighly worse than min-dist clustering method.
        %Turns out, redoing all tree improves SHM correlation. So just redo all
        %trees.

        %Change the root on AncMap
        AncMapCell = calcAncMapCell(Tdata,VDJheader);
        AncMap = AncMapCell{1};
        AncMap(1,2) = size(AncMap,1); %It can be anything, since we'll recalc this.
        AncMap(MinLoc,2) = 0;

        %Update the Root RefSeq.
        Tdata{MinLoc,H.RefSeqLoc} = Tdata{1,H.RefSeqLoc};
        AncMap = sortrows(AncMap,[2 1]);
        Tdata = Tdata(AncMap(:,1),:);

        %Update the refseqs after making the tree
        AncMap = buildTreeLink_B(Tdata,VDJheader);
        Tdata = mapData2AncMap(Tdata,VDJheader,AncMap,'SkipRef');
    
        %Update the child count per parent
        ChildCt = zeros(size(AncMap,1),1);
        for k = 1:size(AncMap,1)
            ChildCt(k) = length(findChild(AncMap,AncMap(k,1)));
        end
        Tdata(:,H.ChildCountLoc) = num2cell(ChildCt);
        
        VDJdata(IdxLoc,:) = Tdata;
    catch
        WarningMsg = sprintf('Warning at %s, sequence group # %d',mfilename,UnqGrpNum(y));
        disp(WarningMsg);
    end
end