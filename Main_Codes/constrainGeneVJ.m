%constrainGeneVJ will correct the V and J lengths based on expected
%location of C and F/W. This is to prevent having a V length shorter than
%the 104C location, and same for 118W for J gene.

function VDJdata = constrainGeneVJ(VDJdata,NewHeader,varargin)
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
    VmapNum = VDJdata{j,FamNumLoc(1)};    
    VallowedDel = Vmap{VmapNum(1),end} - 3; %Subtract 3 since you want to preserve codon of C
    Vdel = VDJdata{j,DelLoc(1)};
    if VallowedDel < 0; VallowedDel = 25; end

    JmapNum = VDJdata{j,FamNumLoc(end)};    
    JallowedDel = Jmap{JmapNum(1),end} - 1; %Subtract 1 since you want to preserve codon of WF
    Jdel = VDJdata{j,DelLoc(end)};
    if JallowedDel < 0; JallowedDel = 25; end

    if Vdel > VallowedDel || Jdel > JallowedDel %Correction Needed
        Seq = VDJdata{j,SeqLoc};
        VMDNJ = cell2mat(VDJdata(j,LengthLoc));
        if Vdel > VallowedDel
            dV = Vdel - VallowedDel; %Force the V length change this much
        else
            dV = 0;
        end
        if Jdel > JallowedDel
            dJ = Jdel - JallowedDel; %Force the J length change this much
        else
            dJ = 0;
        end
        VMDNJ(1) = VMDNJ(1) + dV;
        VMDNJ(end) = VMDNJ(end) + dJ;
        
        if length(Seq) - (VMDNJ(1) + VMDNJ(end)) <= 0
            disp('Cannot readjust C and W loc as this the correction is too large');
            continue
        end
        
        if VMDNJ(1) + VMDNJ(end) >= length(Seq) %Can't work out, since there's no D. 
            disp('Skipping V J correction since this leaves no D behind');
            continue
        end
        
        %Extract the remaining MDN length, format output
        SeqMDN = Seq(VMDNJ(1)+1:end-VMDNJ(end));
        AllowedMiss = ceil(VDJdata{j,SHMLoc(1)} / VMDNJ(1) * length(SeqMDN));
        Dmatch = findGeneMatch(SeqMDN,Dmap,'D',AllowedMiss);
        
        VMDNJ(2:4) = Dmatch{4};
        
        %Update VDJdata with new D match and VMDNJ lengths
        VDJdata(j,LengthLoc) = num2cell(VMDNJ);
        VDJdata(j,DelLoc(1)) = {Vdel-dV};
        VDJdata(j,DelLoc(end)) = {Jdel-dJ};        
        VDJdata(j,FamNumLoc(2)) = Dmatch(1);
        VDJdata(j,FamLoc(2)) = Dmatch(2);
        VDJdata(j,DelLoc(2:3)) = num2cell(Dmatch{3}(1,[1 3]));
        
        %Fill in the details now
        VDJdata(j,:) = buildVDJalignment(VDJdata(j,:),NewHeader,Vmap,Dmap,Jmap); %Alignment Info
        VDJdata(j,:) = makeClassifier(VDJdata(j,:),NewHeader); %Classifier + FormattedSeq
        VDJdata(j,:) = appendMutCt(VDJdata(j,:),NewHeader); %SHM infor on the VMDNJ segments
        VDJdata(j,:) = buildRefSeq(VDJdata(j,:),NewHeader,'germline','single'); %must do singles, since group therapy not done.
        VDJdata(j,:) = findCDR3(VDJdata(j,:),NewHeader); %Get the CDR3 seq and info 
        
        disp('Corrected V J alignment based on C W location');    
    end
end