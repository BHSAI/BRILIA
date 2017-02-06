%fixGeneIndel will correct the V gene inserstion/deletion (indel) error as
%produced by the sequencing platform error. Will not do D and J indel
%corrections since they are generall too short to be certain. 
%
%  VDJdata = fixGeneIndel(VDJdata, VDJheader)
%
%  VDJdata = fixGeneIndel(VDJdata, VDJheader, Vmap, Dmap, Jmap)
%
%  [VDJdata, BadIdx] = fixGeneIndel(VDJdata, VDJheader)

function [VDJdata,varargout] = fixGeneIndel(VDJdata,VDJheader,varargin)
%Extract the VDJ database
if length(varargin) >= 3
    [Vmap,Dmap,Jmap] = deal(varargin{1:3});
else
    [Vmap, Dmap, Jmap] = getCurrentDatabase;
end
H = getHeaderVar(VDJheader);

BadIdx = zeros(size(VDJdata,1),1,'logical');
UpdateIdx = zeros(size(VDJdata,1),1,'logical');
for j = 1:size(VDJdata,1)
    try
        %Correct if deletion is too excessive
        VmapNum = VDJdata{j,H.FamNumLoc(1)};    
        VallowedDel = Vmap{VmapNum(1),end} - 3; %Subtract 3 since you want to preserve codon of C
        if VallowedDel < 0; VallowedDel = 25; end
        Vdel = VDJdata{j,H.DelLoc(1)}; %Correction needed if deletion exceeds allowed deletion

        %Correct if the mismatch count is too excessive
        VmissCt = VDJdata{j,H.VmutLoc};
        Vlen = VDJdata{j,H.LengthLoc(1)};
        Videntity = VmissCt/Vlen; %Normally, VDJ V's are 98% conserved. Hence, 95% identity is questionable. 
        if Vdel > VallowedDel || Videntity >= 0.30 %greater than 30% missed
            Seq = VDJdata{j,H.SeqLoc};

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
            VDJdata{j,H.SeqLoc} = FinalSeq(end-length(Seq)+1:end); %Ensuring same length seq
            VDJdata(j,:) = findVDJmatch(VDJdata(j,:),VDJheader,Vmap,Dmap,Jmap);
            UpdateIdx(j) = 1;
        end
    catch
        ErrorMsg = sprintf('Errored at %s, sequence # %d',mfilename,j);
        disp(ErrorMsg);
        VDJdata{j,H.MiscLoc} = ErrorMsg;
        BadIdx(j) = 1;        
    end
end

%Update those that have changed
VDJdata(UpdateIdx,:) = buildRefSeq(VDJdata(UpdateIdx,:),VDJheader,'germline','first');
VDJdata(UpdateIdx,:) = updateVDJdata(VDJdata(UpdateIdx,:),VDJheader,varargin);

if nargout >=2
    varargout{1} = BadIdx;
end
