%removeDupSeq will look for duplicate sequences, combine them, and remove
%the extra entries. Duplicate sequences can arise after Indel correction or
%after padding/trimming sequences of variable lengths.
%
%  VDJdata = removeDupSeq(VDJdata,NewHeader)

function VDJdata = removeDupSeq(VDJdata,NewHeader)
getHeaderVar;

DelLoc = zeros(size(VDJdata,1),1,'logical');
for j = 1:size(VDJdata,1) - 1
    DupLoc = zeros(size(VDJdata,1),1,'logical');
    Xcount = zeros(size(VDJdata,1),1);

    Seq1 = VDJdata{j,SeqLoc};
    Xloc1 = Seq1 == 'X';
    for k = j+1:size(VDJdata,1)
        Seq2 = VDJdata{k,SeqLoc};
        if DelLoc(k); continue; end
        if length(Seq1) ~= length(Seq2); continue; end
        
        Xloc2 = Seq2 == 'X';
        
        MatchLoc = (Seq1 == Seq2) | Xloc1 | Xloc2;
        if sum(MatchLoc) == length(Seq1)
            DupLoc(k) = 1;
            Xcount(k) = sum(Xloc2);
        end
    end
    
    if max(DupLoc) == 1
        %Find the sequence with the least number of X's
        DupLoc(j) = 1;
        Xcount(j) = sum(Xloc1);
        BestSeqLoc = find(Xcount == min(Xcount(DupLoc)));
        
        %Add up all template counts
        NewTempCt = sum(cell2mat(VDJdata(DupLoc,TemplateLoc)));
        
        %Condense duplicate seq to j, delete all others
        VDJdata{j,SeqLoc} = VDJdata{BestSeqLoc(1),SeqLoc};
        VDJdata{j,TemplateLoc} = NewTempCt;
        
        %Set all duplicate seq past jth seq to empty
        DupLoc(j) = 0;
        DelLoc(DupLoc) = 1;
    end
end
VDJdata(DelLoc,:) = [];


%Check for duplicate nts
% SamSeq = VDJdata(:,SeqLoc);
% [UnqNT,~,UnqNTidx] = unique(SamSeq);% 
% DelIdx = zeros(size(VDJdata,1),1) > 1;
% if max(UnqNTidx) ~= length(SamSeq) %Duplicates found
%     for q = 1:length(UnqNT)
%         SubIdx = find(UnqNTidx == q);
%         if length(SubIdx) > 1
%             TempCt = cell2mat(VDJdata(SubIdx,TemplateLoc));
%             MaxLoc = find(TempCt == max(TempCt));
%             if isempty(MaxLoc)
%                 continue
%             end
%             MaxLoc = MaxLoc(1);
%             MaxIdx = SubIdx(MaxLoc);
%             TotTemp = sum(TempCt);
%             VDJdata{MaxIdx,TemplateLoc} = TotTemp;
%             SubIdx(MaxLoc) = [];
%             DelIdx(SubIdx) = 1;           
%         end
%     end
% end
% 
% VDJdata(DelIdx,:) = [];
        
