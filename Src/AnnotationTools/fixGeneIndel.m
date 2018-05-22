%fixGeneIndel will correct the V gene inserstion/deletion (indel) error as
%produced by the sequencing platform error. Will NOT do D and J indel
%corrections since they are generally too short to be certain. 
%
%  VDJdata = fixGeneIndel(VDJdata, VDJheader, DB)
%
%  [VDJdata, BadIdx] = fixGeneIndel(...)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    VDJheader: main BRILIA header cell
%    DB: Gene database structure (getCurrentDatabase.m)
%
%  OUTPUT
%    VDJdata: modified VDJdata where nt insertion and deletion in the V
%      genes are corrected for. 
%
%  NOTE
%    This code DOES NOT fix indel error that are >= 2 consecutive nts long, 
%      because this could signifiy an alignment error as opposed to an
%      indel error. Also, V and J gene indels are not corrected for.
%
%    In order for indel correction to trigger, it must detect V alignments
%      that either:
%        1) shows a germline V gene deletion that goes beyond the conserved
%           C residue, OR
%        2) has a 20% miss rate in the V alignment caused by indel-induced
%           frame shift error
function VDJdata = fixGeneIndel(VDJdata, Map, DB)
%Extract the VDJ database
M = getMapHeaderVar(DB.MapHeader);

%Begin fixing indel error in V genes
UpdateIdx = zeros(size(VDJdata, 1), 1, 'logical'); %Marks which ones need updates
for k = 1:length(Map.Chain)
    %Determine heavy or light or both IG chain
    if strcmpi(Map.Chain(k), 'H')
        SeqLoc      = Map.hSeq;
        DelLoc      = Map.hDel(1);
        VmutLoc     = Map.hVmut;
        VLengthLoc  = Map.hLength(1);
        MLengthLoc  = Map.hLength(2); 
        GeneNumLoc  = Map.hGeneNum(1);
        GeneNameLoc = Map.hGeneName(1);        
        Vmap = DB.Vmap;
        VkCount = 0; %Set to 0, doesn't do any Vnum shifting.
    elseif strcmpi(Map.Chain(k), 'L')
        SeqLoc      = Map.lSeq;
        DelLoc      = Map.lDel(1);
        VmutLoc     = Map.lVmut;
        VLengthLoc  = Map.lLength(1);
        MLengthLoc  = Map.lLength(2); 
        GeneNumLoc  = Map.lGeneNum(1);
        GeneNameLoc = Map.lGeneName(1);
        Vmap = [DB.Vkmap; DB.Vlmap]; %Easier if it's combined
        VkCount = size(DB.Vkmap, 1); %Need to shift lambda maps by this
    else
        warning('%s: Unknown IG chain %s', mfilename, Chain(k));
        continue
    end
    
    %Pre-determined allowed V ref seq and anchors to avoid broadcasting
    VseqSlice = cell(size(VDJdata, 1), 1); %Vref and Vanchor
    VanchorSlice = cell(size(VDJdata, 1), 1); %Vref and Vanchor
    for j = 1:size(VDJdata, 1)
        Vnum = VDJdata{j, GeneNumLoc};        
        Vname = VDJdata{j, GeneNameLoc};
        if isempty(Vnum); continue; end
        if isempty(Vname); continue; end

        %Determine locus and if lambda, shift Vmap number
        if ~isempty(regexpi(Vname, 'IGLV', 'once')) %Lambda, so shift map num
            Vnum = Vnum + VkCount;
        end
        
        %Construct the sliced Vseq and Vanchor cells
        if ~isempty(Vnum) && Vnum(1) > 0
            VseqSlice{j} = Vmap{Vnum, M.Seq};
            VanchorSlice{j} = Vmap{Vnum, M.Anchor};
        end
    end
    
    %Correct indel
    parfor j = 1:size(VDJdata, 1)
        Tdata = VDJdata(j, :);
        
        %Extract all necessary information
        Seq = Tdata{1, SeqLoc};
        Vdel = Tdata{1, DelLoc};
        Vshm = Tdata{1, VmutLoc};
        Vlen = Tdata{1, VLengthLoc};
        Mlen = Tdata{1, MLengthLoc};
        
        %Make sure all info is available before proceeding
        if isempty(Seq); continue; end
        if isempty(Vdel); continue; end
        if isempty(Vshm); continue; end
        if isempty(Vlen); continue; end
        if isempty(Mlen); continue; end
        if isempty(VanchorSlice{j}); continue; end
        if isempty(VseqSlice{j}); continue; end
        if min([Vshm, Vnum, Vdel, Vlen, Mlen]) < 0; continue; end
        
        %Calc the allowed deletion of V gene 3' end and Valigment miss rate
        VallowedDel = VanchorSlice{j} - 3; %Subtract 3 since you want to preserve codon of C
        if VallowedDel < 0; VallowedDel = 25; end

        %Check if V and Nvd extends past 104C codon
        if Vlen + Mlen < Vlen + Vdel - VallowedDel %Not enough M region for new V
            continue;
        end
        
        %Perform indel check Vdel is too large or Vshm >= 30% of Vlen
        if Vdel > VallowedDel || Vshm/Vlen >= 0.30
            %Reanchor C if Vdel > VallowedDel
            if Vdel > VallowedDel
                AdjustLen = Vdel - VallowedDel;
                Vlen = Vlen + AdjustLen;
                Mlen = Mlen - AdjustLen;
                Vdel = VallowedDel;
            end
            
            %Determine minimum Seq needed to reach RefGene's C codon
            VaddLen = Vdel - VallowedDel;
            VendLoc = Vlen + VaddLen;
            if VendLoc > (Vlen + Mlen) || VendLoc < 1; continue; end %Nonsense VendLoc
            SeqV = Seq(1:VendLoc);

            %Take the full V ref seq up to the same C codon
            SeqVrefFull = VseqSlice{j};
            if VallowedDel-1 >= length(SeqVrefFull); continue; end %Bad Vmap seq
            SeqVref = SeqVrefFull(1:end-VallowedDel);

            %Check for sequencing error insertion/deletion (indel)
            [~, Alignment, StartAt] = swalign(SeqV, SeqVref, 'alphabet', 'nt');
            if max(StartAt) == 0; continue; end %No alignment
            Deletes = regexp(Alignment(1, :), '-'); %Location of deletion error in SeqV
            Inserts = regexp(Alignment(3, :), '-'); %Location of insertion error in SeqV
            if isempty(Deletes) && isempty(Inserts); continue; end %Nothing to change

            %Make sure there're no consecutive del/ins, as these are NOT
            %likely sequencing errors, but rather SHMs or other stuff. 
            ConsecDel = sum(diff(Deletes) == 1);
            ConsecIns = sum(diff(Inserts) == 1);
            if ConsecDel >= 1 || ConsecIns >= 1; continue; end

            %Fix indel in SeqV
            FixedSeq = Alignment(1, :);
            if ~isempty(Deletes) %Fill in the deletions
                FixedSeq(Deletes) = Alignment(3, Deletes);
            end
            if ~isempty(Inserts) %Remove the insertions
                FixedSeq(:, Inserts) = [];
            end

            %If left is trimmed, add trimmed seq to TempSeq (swalign issue)
            if StartAt(1) ~= 1
                FixedSeq = sprintf('%s%s', SeqV(1:StartAt(1)-1), FixedSeq);
            end

            %If right is trimmed, including from swalign & Indel, read to
            %FixedSeq (swalign issue)
            Radd = length(SeqV) - (size(Alignment, 2) + StartAt(1) - 1 - length(Deletes));
            if Radd > 0
                FixedSeq(end+1:end+Radd) = SeqV(end-Radd+1:end);
            end
            
            %Determine if you need to trim 5' end
            if length(FixedSeq) > length(SeqVref) 
                TrimLen = length(FixedSeq) - length(SeqVref);
                FixedSeq(1:TrimLen) = [];
            end

            %Update fields
            NewVlen = Vlen - length(SeqV) + length(FixedSeq);
            if min([NewVlen Mlen]) < 0; continue; end %Invalide correction
            
            FinalSeq = [FixedSeq Seq(VendLoc+1:end)];
            Tdata{1, SeqLoc} = FinalSeq; %New seq
            Tdata{1, VLengthLoc} = NewVlen; %New V len
            Tdata{1, MLengthLoc} = Mlen;
            
            %Mark which one is updated - must redo CDR3, SHM, refSeq
            VDJdata(j, :) = Tdata;
            UpdateIdx(j) = 1;
        end
    end
end

VDJdata(UpdateIdx, :) = buildRefSeq(VDJdata(UpdateIdx, :), Map, DB);
VDJdata(UpdateIdx, :) = updateVDJdata(VDJdata(UpdateIdx, :), Map, DB);
