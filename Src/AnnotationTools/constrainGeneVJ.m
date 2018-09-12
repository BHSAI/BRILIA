%constrainGeneVJ will correct the V and J lengths based on expected
%location of C and F/W. This is to prevent having a V length shorter than
%the 104C location, and same for 118W for J gene. If there is an adjustment
%that is required, it will realign the D segment to get a new Nvd-D-Ndj
%segment.
%
%  VDJdata = constrainGeneVJ(VDJdata, VDJheader, DB)
%
%  [VDJdata, BadIdx] = constrainGeneVJ(...)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    VDJheader: main BRILIA header cell
%    DB: Gene database structure (getCurrentDatabase.m)
%
%  OUTPUT
%    VDJdata: modified VDJdata where V and J genes are set to cover the
%      CDR3 region

function VDJdata = constrainGeneVJ(VDJdata, Map, DB)
if isempty(VDJdata); return; end
%Determine chain and extract key locations
M = getMapHeaderVar(DB.MapHeader);

%Check each entry to ensure V and J covers CDR3
UpdateIdx = zeros(size(VDJdata, 1), 1, 'logical');
for k = 1:length(Map.Chain) 
    if strcmpi(Map.Chain(k), 'H')
        DelLoc  = Map.hDel;
        GeneNumLoc  = Map.hGeneNum; 
        GeneNameLoc  = Map.hGeneName; 
        SeqLoc  = Map.hSeq;
        LengthLoc  = Map.hLength;
        Vmap = DB.Vmap;
        Jmap = DB.Jmap;
        Dmap = DB.Dmap;
    else
        DelLoc  = Map.lDel;
        GeneNumLoc  = Map.lGeneNum; 
        GeneNameLoc  = Map.lGeneName; 
        SeqLoc  = Map.lSeq;
        LengthLoc  = Map.lLength;
        Vmap = [DB.Vkmap; DB.Vlmap];
        Jmap = [DB.Jkmap; DB.Jlmap];
        VkCount = size(DB.Vkmap, 1);
        JkCount = size(DB.Jkmap, 1);
    end
    
    for j = 1:size(VDJdata, 1)
        %Extract necessary info
        VmapNum = VDJdata{j, GeneNumLoc(1)};
        JmapNum = VDJdata{j, GeneNumLoc(end)};
        Vname = VDJdata{j, GeneNameLoc(1)};
        Jname = VDJdata{j, GeneNameLoc(end)};
        Vdel = VDJdata{j, DelLoc(1)};
        Jdel = VDJdata{j, DelLoc(end)};
        Seq = VDJdata{j, SeqLoc};
        SegLen = cell2mat(VDJdata(j, LengthLoc));
        
        %Make sure all necessary info is there
        if isempty(VmapNum) 
            continue; 
        else
            VmapNum = VmapNum(1);
        end
        if isempty(JmapNum) 
            continue; 
        else
            JmapNum = JmapNum(1);
        end
        if isempty(Vname); continue; end
        if isempty(Jname); continue; end
        if isempty(Vdel); continue; end
        if isempty(Jdel); continue; end
        if isempty(Seq); continue; end
        if isempty(SegLen); continue; end
        if min([Vdel(:); Jdel(:); SegLen(:)]) < 0; continue; end
        
        %If it's a light chain lambda, need to shift Map num
        if ~isempty(regexpi(Vname, 'IGL[VDJ]', 'once'))
            VmapNum = VmapNum + VkCount;
            JmapNum = JmapNum + JkCount;
        end
        
        %Calc max Vref del allowed to retain 104C
        VallowedDel = Vmap{VmapNum, M.Anchor} - 3; %Subtract 3 since you want to preserve codon of C
        if VallowedDel < 0; VallowedDel = 25; end
        
        %Calc max Jref del allowed to retain 118W
        JallowedDel = Jmap{JmapNum, M.Anchor} - 1; %Subtract 1 since you want to preserve codon of WF
        if JallowedDel < 0; JallowedDel = 25; end
        
        %Calc Vdel and Jdel adjustments
        dV = 0;
        if Vdel > VallowedDel
            dV = Vdel - VallowedDel; %Force the V length change this much
        end
        dJ = 0;
        if Jdel > JallowedDel
            dJ = Jdel - JallowedDel; %Force the J length change this much
        end
        
        %Update VDJdata if you need to adjust ref deletions
        if dJ > 0 || dV > 0 %Correction Needed
            SegLen(1) = SegLen(1) + dV;
            SegLen(end) = SegLen(end) + dJ;

            if SegLen(1) + SegLen(end) >= length(Seq)
                %disp('Skipping V J correction since this leaves no D behind');
                continue
            end
            
            %Recalculate new V and J dels
            NewVdel = Vdel - dV;
            NewJdel = Jdel - dJ;
    
            %Redo the D alignment for heavy chain
            if Map.Chain(k) == 'H'
                SeqMDN = Seq(SegLen(1)+1:end-SegLen(end));
                Dmatch = findGeneMatch(SeqMDN, Dmap, 'D', 0);
                SegLen(2:4) = Dmatch{1, 4};
                VDJdata(j, [GeneNumLoc(2) GeneNameLoc(2)]) = Dmatch(1, 1:2);
                VDJdata(j, DelLoc) = {NewVdel Dmatch{1, 3}(1) Dmatch{1, 3}(3) NewJdel};
                VDJdata(j, LengthLoc) = num2cell(SegLen);                
            else
                SegLen(2) = length(Seq) - SegLen(1) - SegLen(end);
                VDJdata(j, DelLoc) = {NewVdel NewJdel};
                VDJdata(j, LengthLoc) = num2cell(SegLen);                
            end

            %Mark entries that have been updated
            UpdateIdx(j) = 1;
        end
    end
end

%Update those that have changed
if max(UpdateIdx) > 0 
    VDJdata(UpdateIdx, :) = buildRefSeq(VDJdata(UpdateIdx, :), Map, DB, Map.Chain, 'germline', 'first');
    VDJdata(UpdateIdx, :) = updateVDJdata(VDJdata(UpdateIdx, :), Map, DB);
    showStatus(sprintf('  Corrected %d V or J ends', sum(UpdateIdx)));
end
