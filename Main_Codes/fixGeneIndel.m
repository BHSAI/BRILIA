%Correct the V sequence error deletion/insertion locations
function VDJdata = fixGeneIndel(VDJdata,NewHeader,varargin)
%Extract the VDJ database
if length(varargin) == 3
    Vmap = varargin{1};
    Dmap = varargin{2};
    Jmap = varargin{3};
else
    [Vmap, Dmap, Jmap] = getCurrentDatabase;
end

getHeaderVar;

for j = 1:size(VDJdata,1)
    %Correct if deletion is too excessive
    VmapNum = VDJdata{j,FamNumLoc(1)};    
    VallowedDel = Vmap{VmapNum(1),end} - 3; %Subtract 3 since you want to preserve codon of C
    if VallowedDel < 0; VallowedDel = 25; end
    Vdel = VDJdata{j,DelLoc(1)}; %Correction needed if deletion exceeds allowed deletion
    
    %Correct if the mismatch count is too excessive
    VmissCt = VDJdata{j,SHMLoc(1)};
    Vlen = VDJdata{j,LengthLoc(1)};
    Videntity = VmissCt/Vlen; %Normally, VDJ V's are 98% conserved. Hence, 95% identity is questionable. 
    if Vdel > VallowedDel || Videntity >= 0.10 %greater than 10% missed
        Seq = VDJdata{j,SeqLoc};

        %Determine minimum Seq needed to reach RefGene's C codon
        VaddLen = Vdel - VallowedDel;
        VendLoc = Vlen(1) + VaddLen;
        SeqV = Seq(1:VendLoc);

        %Take the full ref seq
        SeqVrefFull = Vmap{VmapNum(1),1};
        SeqVref = SeqVrefFull(1:end - VallowedDel);        
        if isempty(SeqV); continue; end %Broken sequence, so skip it
        [~, Alignment, StartAt] = swalign(SeqV,SeqVref,'alphabet','nt');
        
        %Check for sequencing error insertion/deletion (indel), as that could happen        
        Deletes = regexp(Alignment(1,:),'-'); %Location of deletion error
        Inserts = regexp(Alignment(3,:),'-'); %Location of insertion error
        
        if isempty(Deletes) && isempty(Inserts); continue; end

        %Make sure there's no consectuvie del/ins, as this is unlikely to
        %be caused by sequencing error - they normablly do single indels
        %here and there.
        ConsecDel = sum(diff(Deletes) == 1);
        ConsecIns = sum(diff(Inserts) == 1);
        if ConsecDel > 1 || ConsecIns > 1; continue; end

        TempSeq = Alignment(1,:);
        
        %Fill in the deletions
        if ~isempty(Deletes)
            TempSeq(Deletes) = Alignment(3,Deletes);
        end
        
        %Remove the insertions
        if ~isempty(Inserts)
            TempSeq(:,Inserts) = [];
        end
        
        %If left is trimmed, read to TempSeq
        if StartAt(1) ~= 1
            TempSeq = [SeqV(1:StartAt(1)-1) TempSeq];
        end
        
        %If right is trimmed, including from swalign & Indel, readd to TempSeq
        Radd = length(SeqV) - (size(Alignment,2) + StartAt(1) - 1 - length(Deletes));
        if Radd > 0
            TempSeq(end+1:end+Radd) = SeqV(end-Radd+1:end);
        end
        
        %Determine how many NTs' to pull.
        PullLength = length(SeqV) - length(TempSeq);
        if PullLength > 0
            S1 = length(SeqVref) - length(TempSeq) - PullLength + 1;
            S2 = length(SeqVref) - length(TempSeq);
            TempSeq = [SeqVref(S1:S2) TempSeq]; %This is done in case there are so many insertions
        end
        
        %Preserve Seq Length for all seq, trim ends
        FinalSeq = [TempSeq Seq(VendLoc+1:end)];        
        VDJdata{j,SeqLoc} = FinalSeq(end-length(Seq)+1:end); %Ensuring same length seq
        VDJdata(j,:) = findVDJmatch(VDJdata(j,:),NewHeader,Vmap,Dmap,Jmap);
    end
end