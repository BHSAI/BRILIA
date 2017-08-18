%buildVDJalignment will return the VDJ alignment results for each sequence.
%An "-" is a match with the germline gene, a lower case are deleted nts
%from the germline gene, and a capital letter is a mismatch with the
%germline gene, and an "|" marks where the germline gene was deleted.
%
%  VDJdata = buildVDJalignment(VDJdata, VDJheader)
%
%  VDJdata = buildVDJalignment(VDJdata, VDJheader, DB)

function [VDJdata, VDJheader] = buildVDJalignment(VDJdata,VDJheader,DB)
[H, L, Chain] = getAllHeaderVar(VDJheader);
   
if ~isempty(strfind(Chain, 'H'))
    %Build the heavy chain alignment
    HeavyAlign = cell(size(VDJdata,1),3);
    for j = 1:size(VDJdata,1)
        %Get all the necessary data
        Seq = VDJdata{j,H.SeqLoc};
        VMDNJ = cell2mat(VDJdata(j,H.LengthLoc));

        %Make sure all entries are valid
        if isempty(Seq);  continue; end
        if VMDNJ(1) <= 0; continue; end
        if VMDNJ(2) <  0; continue; end
        if VMDNJ(3) <= 0; continue; end
        if VMDNJ(4) <  0; continue; end
        if VMDNJ(5) <= 0; continue; end

        %----------------------------------------------------------------------
        %Setting V alignment
        Vnum = VDJdata{j,H.GeneNumLoc(1)};
        if ~isempty(Vnum)
            Vnum = Vnum(1);
            VsamSeq = Seq(1:VMDNJ(1));
            VrefSeq = DB.Vmap{Vnum,1};
            VntDel = VDJdata{j,H.DelLoc(1)};

            %Determine if there are extra nts left for padding
            ExtraLeft = VMDNJ(1)+ VntDel - length(VrefSeq);
            if ExtraLeft > 0
                VrefSeq = sprintf('%s%s',repmat('X',ExtraLeft,1),VrefSeq);
            elseif ExtraLeft < 0 %Delete some
                VrefSeq(1:abs(ExtraLeft)) = [];
            end

            %Get the deleted nts
            if VntDel > 0
                Vref3 = VrefSeq(end-VntDel+1:end);
                VrefSeq(end-VntDel+1:end) = [];
            else
                Vref3 = '';
            end

            %Assembe the V alignment string and save
            Valign = VrefSeq;
            Valign(VsamSeq == VrefSeq | VsamSeq == 'X') = '.';
            HeavyAlign{j,1} = sprintf('%s|%s',Valign,lower(Vref3));
        end

        %----------------------------------------------------------------------
        %Setting D alignment
        Dnum = VDJdata{j,H.GeneNumLoc(2)};
        if ~isempty(Dnum)
            Dnum = Dnum(1);
            DsamSeq = Seq(sum(VMDNJ(1:2))+1:sum(VMDNJ(1:3)));

            %Get the deleted nts
            DntDel5 = VDJdata{j,H.DelLoc(2)};
            DntDel3 = VDJdata{j,H.DelLoc(3)};        
            DrefSeq = DB.Dmap{Dnum,1}(1+DntDel5:end-DntDel3);       
            if DntDel5 > 0
                Dref5 = DB.Dmap{Dnum,1}(1:DntDel5);
            else
                Dref5 = '';
            end
            if DntDel3 > 0
                Dref3 = DB.Dmap{Dnum,1}(end-DntDel3+1:end);
            else
                Dref3 = '';
            end

            %Assemble the alignment and save
            Dalign = DrefSeq;
            Dalign(DsamSeq==DrefSeq | DsamSeq == 'X') = '.';
            HeavyAlign{j,2} = sprintf('%s|%s|%s',lower(Dref5),Dalign,lower(Dref3));
        end

        %----------------------------------------------------------------------
        %Setting J alignment
        Jnum = VDJdata{j,H.GeneNumLoc(3)};
        if ~isempty(Jnum)
            Jnum = Jnum(1);
            JsamSeq = Seq(end-VMDNJ(5)+1:end);
            JrefSeq = DB.Jmap{Jnum,1};
            JntDel = VDJdata{j,H.DelLoc(4)};

            %Determine if there are extra nts left for padding
            ExtraRight = VMDNJ(5) + JntDel - length(JrefSeq);
            if ExtraRight > 0
                JrefSeq = sprintf('%s%s',JrefSeq,repmat('X',ExtraRight,1));
            elseif ExtraRight < 0
                JrefSeq(end-abs(ExtraRight)+1:end)= [];
            end

            %Get the deleted nts
            if JntDel > 0
                Jref5 = JrefSeq(1:JntDel);
                JrefSeq(1:JntDel) = [];
            else
                Jref5 = '';
            end

            %Assemble the alignment and save
            Jalign = JrefSeq;
            Jalign(JsamSeq == JrefSeq | JsamSeq == 'X') = '.';
            HeavyAlign{j,3} = sprintf('%s|%s',lower(Jref5),Jalign);
        end
    end
    
    %Append or replace alignment in VDJdata
    HeavyAlignLoc = findCell(VDJheader,{'H-V_Align','H-D_Align','H-J_Align'});
    if max(HeavyAlignLoc == 0) == 1 %If anything is 0
        if sum(HeavyAlignLoc > 0) ~= 3 %Something went wrong and not all alignments were done. delete and redo.
            VDJdata(:,HeavyAlignLoc(HeavyAlignLoc>0)) = [];
        end
        VDJdata = [VDJdata HeavyAlign];
        VDJheader = [VDJheader 'H-V_Align' 'H-D_Align' 'H-J_Align'];
    else %Already align, just replace
        VDJdata(:,HeavyAlignLoc) = HeavyAlign;
    end
end

if ~isempty(strfind(Chain, 'L'))
    %Build the heavy chain alignment
    LightAlign = cell(size(VDJdata,1),2);
    for j = 1:size(VDJdata,1)
        %Get all the necessary data
        Seq = VDJdata{j,L.SeqLoc};
        VNJ = cell2mat(VDJdata(j,L.LengthLoc));

        %Make sure all entries are valid
        if isempty(Seq);  continue; end
        if VNJ(1) <= 0; continue; end
        if VNJ(2) <  0; continue; end
        if VNJ(3) <= 0; continue; end

        %----------------------------------------------------------------------
        %Setting V alignment
        Vnum = VDJdata{j,L.GeneNumLoc(1)};
        Vname = VDJdata{j,L.GeneNameLoc(1)};
        if ~isempty(Vnum)
            Vnum = Vnum(1);
            VsamSeq = Seq(1:VNJ(1));
            VntDel = VDJdata{j,L.DelLoc(1)};
            
            %Choose the kappa/lampda DB
            if ~isempty(regexpi(Vname,'IGK')) 
                VrefSeq = DB.Vkmap{Vnum,1};
            else
                VrefSeq = DB.Vlmap{Vnum,1};
            end

            %Determine if there are extra nts left for padding
            ExtraLeft = VNJ(1)+ VntDel - length(VrefSeq);
            if ExtraLeft > 0
                VrefSeq = sprintf('%s%s',repmat('X',ExtraLeft,1),VrefSeq);
            elseif ExtraLeft < 0 %Delete some
                VrefSeq(1:abs(ExtraLeft)) = [];
            end

            %Get the deleted nts
            if VntDel > 0
                Vref3 = VrefSeq(end-VntDel+1:end);
                VrefSeq(end-VntDel+1:end) = [];
            else
                Vref3 = '';
            end

            %Assembe the V alignment string and save
            Valign = VrefSeq;
            Valign(VsamSeq == VrefSeq | VsamSeq == 'X') = '.';
            LightAlign{j,1} = sprintf('%s|%s',Valign,lower(Vref3));
        end

        %----------------------------------------------------------------------
        %Setting J alignment
        Jnum = VDJdata{j,L.GeneNumLoc(2)};
        Jname = VDJdata{j,L.GeneNameLoc(2)};
        if ~isempty(Jnum)
            Jnum = Jnum(1);
            JsamSeq = Seq(end-VNJ(3)+1:end);
            JntDel = VDJdata{j,L.DelLoc(2)};
            
            %Choose the kappa/lampda DB
            if ~isempty(regexpi(Jname,'IGK')) 
                JrefSeq = DB.Jkmap{Jnum,1};
            else
                JrefSeq = DB.Jlmap{Jnum,1};
            end
            
            %Determine if there are extra nts left for padding
            ExtraRight = VNJ(3) + JntDel - length(JrefSeq);
            if ExtraRight > 0
                JrefSeq = sprintf('%s%s',JrefSeq,repmat('X',ExtraRight,1));
            elseif ExtraRight < 0
                JrefSeq(end-abs(ExtraRight)+1:end)= [];
            end

            %Get the deleted nts
            if JntDel > 0
                Jref5 = JrefSeq(1:JntDel);
                JrefSeq(1:JntDel) = [];
            else
                Jref5 = '';
            end

            %Assemble the alignment and save
            Jalign = JrefSeq;
            Jalign(JsamSeq == JrefSeq | JsamSeq == 'X') = '.';
            LightAlign{j,2} = sprintf('%s|%s',lower(Jref5),Jalign);
        end
    end
    
    LightAlignLoc = findCell(VDJheader,{'L-V_Align','L-J_Align'});
    if max(LightAlignLoc == 0) == 1 %If anything is 0
        if sum(LightAlignLoc > 0) ~= 2 %Something went wrong and not all alignments were done. delete and redo.
            VDJdata(:,LightAlignLoc(LightAlignLoc>0)) = [];
        end
        VDJdata = [VDJdata LightAlign];
        VDJheader = [VDJheader 'L-V_Align' 'L-J_Align'];
    else %Already align, just replace
        VDJdata(:,LightAlignLoc) = LightAlign;
    end
end

