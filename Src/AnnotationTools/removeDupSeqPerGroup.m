function VDJdata = removeDupSeqPerGroup(VDJdata, Map)
%Look for duplicate entries
DelIdx = zeros(size(VDJdata, 1), 1, 'logical'); %Keeps track of duplicate entries that will be removed
for j = 1:size(VDJdata, 1) - 1
    if DelIdx(j); continue; end %Already deleted, skip
    DupLoc = zeros(size(VDJdata, 1), 1, 'logical'); %Tracks location of duplicate sequence
    Xcount = zeros(size(VDJdata, 1), 1); %Tracks how many X's there are
            
    %Add on constant regions to prevent collapsing seq with different C.
    if contains(Map.Chain, 'H') && ~isempty(VDJdata{j, Map.hOverSeq3})
        ConstSeq1 = VDJdata{j, Map.hOverSeq3};
    else 
        ConstSeq1 = '';
    end

    %Get Seq1 info and location of X
    if strcmpi(Map.Chain, 'HL')
        Seq1 = sprintf('%s%s%s', VDJdata{j, Map.hSeq}, VDJdata{j, Map.lSeq}, ConstSeq1);
    elseif strcmpi(Map.Chain, 'H')
        Seq1 = sprintf('%s%s', VDJdata{j, Map.hSeq}, ConstSeq1);
    else
        Seq1 = VDJdata{j, Map.lSeq};
    end
       
    if isempty(Seq1); continue; end
    Xloc1 = Seq1 == 'N';

    %Look for other sequence that are same as Seq1
    for q = j+1:size(VDJdata, 1)
        if DelIdx(q); continue; end %Already deleted, skip

        %Determine if there are any possible constant region sequences
        if contains(Map.Chain, 'H') && ~isempty(VDJdata{q, Map.hOverSeq3})
            ConstSeq2 = VDJdata{q, Map.hOverSeq3};
        else 
            ConstSeq2 = '';
        end
    
        %Get Seq2 info and location of X
        if strcmpi(Map.Chain, 'HL')
            Seq2 = sprintf('%s%s%s', VDJdata{q, Map.hSeq}, VDJdata{q, Map.lSeq}, ConstSeq2);
        elseif strcmpi(Map.Chain, 'H')
            Seq2 = sprintf('%s%s', VDJdata{q, Map.hSeq}, ConstSeq2);
        else
            Seq2 = VDJdata{q, Map.lSeq};
        end
        if length(Seq1) ~= length(Seq2); continue; end %Can't be same
        Xloc2 = Seq2 == 'N';

        %See if they are the same
        MatchIdx = (Seq1 == Seq2) | Xloc1 | Xloc2;
        if sum(MatchIdx) == length(Seq1)
            DupLoc(q) = 1;
            Xcount(q) = sum(Xloc2);
        end
    end

    if max(DupLoc)
        %Find the sequence with the least number of X's
        DupLoc(j) = 1; %Fill in 1st entry that was skipped before
        Xcount(j) = sum(Xloc1); %Fill in 1st entry
        BestSeqLoc = find((Xcount == min(Xcount(DupLoc))) & DupLoc); %get the one with least X's
        BestSeqLoc = BestSeqLoc(1); %Get the first best one

        %Add up all template counts
        NewTempCt = sum(cell2mat(VDJdata(DupLoc, Map.Template)));

        %Condense duplicate seq to j, delete all others
        VDJdata(j, :) = VDJdata(BestSeqLoc, :);
        VDJdata{j, Map.Template} = NewTempCt;

        %Set all duplicate seq past jth seq to empty
        DupLoc(j) = 0;
        DelIdx(DupLoc) = 1;
    end
end
VDJdata(DelIdx, :) = [];