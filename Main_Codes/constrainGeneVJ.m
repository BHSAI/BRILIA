%constrainGeneVJ will correct the V and J lengths based on expected
%location of C and F/W. This is to prevent having a V length shorter than
%the 104C location, and same for 118W for J gene. If there is an adjustment
%that is required, it will realign the D segment to get a new Nvd-D-Ndj
%segment.
%
%  VDJdata = constrainGeneVJ(VDJdata,NewHeader)
%
%  VDJdata = constrainGeneVJ(VDJdata,NewHeader,Vmap,Dmap,Jmap)

function [VDJdata,varargout] = constrainGeneVJ(VDJdata,NewHeader,varargin)
%Extract the VDJ database
if length(varargin) >= 3
    [Vmap,Dmap,Jmap] = deal(varargin{1:3});
else
    [Vmap,Dmap,Jmap] = getCurrentDatabase;
end
getHeaderVar;

BadIdx = zeros(size(VDJdata,1),'logical');
UpdateIdx = zeros(size(VDJdata,1),'logical');
for j = 1:size(VDJdata,1)
    try
        %Calc max Vref del allowed to retain 104C
        VmapNum = VDJdata{j,FamNumLoc(1)}(1);    
        VallowedDel = Vmap{VmapNum,10} - 3; %Subtract 3 since you want to preserve codon of C
        Vdel = VDJdata{j,DelLoc(1)};
        if VallowedDel < 0; VallowedDel = 25; end
        
        %Calc max Jref del allowed to retain 118W
        JmapNum = VDJdata{j,FamNumLoc(end)}(1);    
        JallowedDel = Jmap{JmapNum,10} - 1; %Subtract 1 since you want to preserve codon of WF
        Jdel = VDJdata{j,DelLoc(end)};
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
            Seq = VDJdata{j,SeqLoc};
            VMDNJ = cell2mat(VDJdata(j,LengthLoc));
            VMDNJ(1) = VMDNJ(1) + dV;
            VMDNJ(end) = VMDNJ(end) + dJ;            

            if VMDNJ(1) + VMDNJ(end) >= length(Seq)
                %disp('Skipping V J correction since this leaves no D behind');
                continue
            end

            %Redo the D alignment
            SeqMDN = Seq(VMDNJ(1)+1:end-VMDNJ(end));
            Dmatch = findGeneMatch(SeqMDN,Dmap,'D',0);
            VMDNJ(2:4) = Dmatch{1,4};

            %Update VDJdata with new D match and VMDNJ lengths
            VDJdata(j,LengthLoc) = num2cell(VMDNJ);
            VDJdata(j,DelLoc) = {(Vdel-dV) Dmatch{1,3}(1) Dmatch{1,3}(3) (Jdel-dJ)};
            VDJdata(j,[FamNumLoc(2) FamLoc(2)]) = Dmatch(1,1:2);
            UpdateIdx(j) = 1;
        end
    catch
        ErrorMsg = sprintf('Errored at %s, sequence # %d',mfilename,j);
        disp(ErrorMsg);
        VDJdata{j,MiscLoc} = ErrorMsg;
        BadIdx(j) = 1;        
    end
end

%Update those that have changed
VDJdata(UpdateIdx,:) = buildRefSeq(VDJdata(UpdateIdx,:),NewHeader,'germline','first');
VDJdata(UpdateIdx,:) = updateVDJdata(VDJdata(UpdateIdx,:),NewHeader,varargin);


if nargout >=2
    varargout{1} = BadIdx;
end