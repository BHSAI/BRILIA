%buildVDJalignment will take VDJ data, and return the alignment
%information.

function VDJdata = buildVDJalignment(VDJdata,NewHeader,varargin)
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
    Seq = VDJdata{j,SeqLoc};
    VMDNJ = cell2mat(VDJdata(j,LengthLoc));
    
    %----------------------------------------------------------------------
    %Setting V alignment
    Vnum = VDJdata{j,FamNumLoc(1)};
    if ~isempty(Vnum)
        VsamSeq = Seq(1:VMDNJ(1));
        
        %Extract the germline seq
        Vnum = Vnum(1,1); %Get only the first one
        VntDel = VDJdata{j,DelLoc(1)};
        VrefSeq = Vmap{Vnum,1}(end-VntDel-VMDNJ(1)+1:end-VntDel);        
        if VntDel > 0
            Vref3 = Vmap{Vnum,1}(end-VntDel+1:end);
        else
            Vref3 = '';
        end
        
        Valign = VrefSeq;
        Valign(VsamSeq == VrefSeq) = '-';
        Valign = [Valign '|' lower(Vref3)];
        VDJdata{j,AlignLoc(1)} = Valign;
    end
    
    %----------------------------------------------------------------------
    %Setting D alignment
    Dnum = VDJdata{j,FamNumLoc(2)};
    if ~isempty(Dnum)
        DsamSeq = Seq(sum(VMDNJ(1:2))+1:sum(VMDNJ(1:3)));
        
        %Extract the germline seq
        Dnum = Dnum(1,1); %Get only the first one
        DntDel5 = VDJdata{j,DelLoc(2)};
        DntDel3 = VDJdata{j,DelLoc(3)};        
        DrefSeq = Dmap{Dnum,1}(1+DntDel5:end-DntDel3);       
        if DntDel5 > 0
            Dref5 = Dmap{Dnum,1}(1:DntDel5);
        else
            Dref5 = '';
        end
        if DntDel3 > 0
            Dref3 = Dmap{Dnum,1}(end-DntDel3+1:end);
        else
            Dref3 = '';
        end
        
        Dalign = DrefSeq;
        Dalign(DrefSeq==DsamSeq) = '-';
        Dalign = [lower(Dref5) '|' Dalign '|' lower(Dref3)];
        VDJdata{j,AlignLoc(2)} = Dalign;
    end
    
    %----------------------------------------------------------------------
    %Setting J alignment
    Jnum = VDJdata{j,FamNumLoc(3)};
    if ~isempty(Jnum)
        JsamSeq = Seq(end-VMDNJ(5)+1:end);
        
        %Extract the germline seq
        Jnum = Jnum(1,1); %Get only the first one
        JntDel = VDJdata{j,DelLoc(4)};
        JrefSeq = Jmap{Jnum,1}(1+JntDel:1+JntDel+VMDNJ(5)-1);        
        if JntDel > 0
            Jref5 = Jmap{Jnum,1}(1:JntDel);
        else
            Jref5 = '';
        end
        
        Jalign = JrefSeq;
        Jalign(JsamSeq == JrefSeq) = '-';
        Jalign = [lower(Jref5) '|' Jalign];
        VDJdata{j,AlignLoc(3)} = Jalign;
    end
end