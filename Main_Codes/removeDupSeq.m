%removeDupSeq will look for duplicate sequences, combine them, and remove
%the extra entries. Duplicate sequences can arise after Indel correction or
%after padding/trimming sequences of variable lengths.
%
%  VDJdata = removeDupSeq(VDJdata,VDJheader)

function VDJdata = removeDupSeq(VDJdata,VDJheader)
H = getHeaderVar(VDJheader);

%Look for duplicate entries
H.DelLoc = zeros(size(VDJdata,1),1,'logical'); %Keeps track of duplicate entries that will be removed
for j = 1:size(VDJdata,1) - 1
    if H.DelLoc(j); continue; end %Already deleted, skip
    DupLoc = zeros(size(VDJdata,1),1,'logical'); %Tracks location of duplicate sequence
    Xcount = zeros(size(VDJdata,1),1); %Tracks how many X's there are

    %Get Seq1 info and location of X
    Seq1 = VDJdata{j,H.SeqLoc};
    Xloc1 = Seq1 == 'X';

    %Look for other sequence that are same as Seq1
    for k = j+1:size(VDJdata,1)
        if H.DelLoc(k); continue; end %Already deleted, skip

        %Get Seq2 info and location of X
        Seq2 = VDJdata{k,H.SeqLoc};
        if length(Seq1) ~= length(Seq2); continue; end %Can't be same
        Xloc2 = Seq2 == 'X';
        
        %See if they are the same
        H.MatchLoc = (Seq1 == Seq2) | Xloc1 | Xloc2;
        if sum(H.MatchLoc) == length(Seq1)
            DupLoc(k) = 1;
            Xcount(k) = sum(Xloc2);
        end
    end
    
    if max(DupLoc) == 1
        %Find the sequence with the least number of X's
        DupLoc(j) = 1;
        Xcount(j) = sum(Xloc1);
        BestSeqLoc = find((Xcount == min(Xcount(DupLoc))) & DupLoc);
        BestSeqLoc = BestSeqLoc(1); 
        
        %Add up all template counts
        NewTempCt = sum(cell2mat(VDJdata(DupLoc,H.TemplateLoc)));
        
        %Condense duplicate seq to j, delete all others
        VDJdata(j,:) = VDJdata(BestSeqLoc,:);
        VDJdata{j,H.TemplateLoc} = NewTempCt;
        
        %Set all duplicate seq past jth seq to empty
        DupLoc(j) = 0;
        H.DelLoc(DupLoc) = 1;
    end
end
VDJdata(H.DelLoc,:) = [];
disp(['Removed ' num2str(sum(H.DelLoc)) ' duplicate sequences.']);
