%buildVDJalignment will return the VDJ alignment results for each sequence.
%An "-" is a match with the germline gene, a lower case are deleted nts
%from the germline gene, and a capital letter is a mismatch with the
%germline gene, and an "|" marks where the germline gene was deleted.
%
%  VDJdata = buildVDJalignment(VDJdata, VDJheader)
%
%  VDJdata = buildVDJalignment(VDJdata, VDJheader, DB)

function VDJdata = buildVDJalignment(VDJdata,Map,DB)
if any(contains(Map.Chain, 'h', 'ignorecase', true)) 
    %Build the heavy chain alignment
    HeavyAlign = cell(size(VDJdata,1),3);
    for j = 1:size(VDJdata,1)
        Seq = VDJdata{j,Map.hSeq};
        VMDNJ = cell2mat(VDJdata(j,Map.hLength));
        if isempty(Seq);  continue; end
        if any(VMDNJ < [1 0 1 0 1]); continue; end

        %----------------------------------------------------------------------
        %Setting V alignment
        Vnum = VDJdata{j,Map.hGeneNum(1)};
        if ~isempty(Vnum)
            Vnum = Vnum(1);
            VsamSeq = Seq(1:VMDNJ(1));
            VrefSeq = DB.Vmap{Vnum,1};
            Vdel = VDJdata{j,Map.hDel(1)};

            %Determine if there are extra nts left for padding
            ExtraLeft = VMDNJ(1)+ Vdel - length(VrefSeq);
            if ExtraLeft > 0
                VrefSeq = sprintf('%s%s',repmat('N',ExtraLeft,1),VrefSeq);
            elseif ExtraLeft < 0 %Delete some
                VrefSeq(1:abs(ExtraLeft)) = [];
            end

            %Get the deleted nts
            if Vdel > 0
                Vref3 = VrefSeq(end-Vdel+1:end);
                VrefSeq(end-Vdel+1:end) = [];
            else
                Vref3 = '';
            end

            %Assemble the V alignment string and save
            Valign = VrefSeq;
            Valign(VsamSeq == VrefSeq | VsamSeq == 'N') = '.';
            HeavyAlign{j,1} = sprintf('%s|%s',Valign,lower(Vref3));
        end

        %----------------------------------------------------------------------
        %Setting D alignment
        Dnum = VDJdata{j,Map.hGeneNum(2)};
        if ~isempty(Dnum)
            Dnum = Dnum(1);
            DsamSeq = Seq(sum(VMDNJ(1:2))+1:sum(VMDNJ(1:3)));

            %Get the deleted nts
            DntDel5 = VDJdata{j,Map.hDel(2)};
            DntDel3 = VDJdata{j,Map.hDel(3)};        
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
            Dalign(DsamSeq==DrefSeq | DsamSeq == 'N') = '.';
            HeavyAlign{j,2} = sprintf('%s|%s|%s',lower(Dref5),Dalign,lower(Dref3));
        end

        %----------------------------------------------------------------------
        %Setting J alignment
        Jnum = VDJdata{j,Map.hGeneNum(3)};
        if ~isempty(Jnum)
            Jnum = Jnum(1);
            JsamSeq = Seq(end-VMDNJ(5)+1:end);
            JrefSeq = DB.Jmap{Jnum,1};
            JntDel = VDJdata{j,Map.hDel(4)};

            %Determine if there are extra nts left for padding
            ExtraRight = VMDNJ(5) + JntDel - length(JrefSeq);
            if ExtraRight > 0
                JrefSeq = sprintf('%s%s',JrefSeq,repmat('N',ExtraRight,1));
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
            Jalign(JsamSeq == JrefSeq | JsamSeq == 'N') = '.';
            HeavyAlign{j,3} = sprintf('%s|%s',lower(Jref5),Jalign);
        end
    end
    
    %Append or replace alignment in VDJdata
    HeavyAlignLoc = [Map.hValign Map.hDalign Map.hJalign];
    VDJdata(:,HeavyAlignLoc) = HeavyAlign;
end

if any(contains(Map.Chain, 'l', 'ignorecase', true))
    %Build the heavy chain alignment
    LightAlign = cell(size(VDJdata,1),2);
    for j = 1:size(VDJdata,1)
        Seq = VDJdata{j,Map.lSeq};
        VNJ = cell2mat(VDJdata(j,Map.lLength));
        if isempty(Seq);  continue; end
        if any(VNJ < [1 0 1]); continue; end

        %----------------------------------------------------------------------
        %Setting V alignment
        Vnum = VDJdata{j,Map.lGeneNum(1)};
        Vname = VDJdata{j,Map.lGeneName(1)};
        if ~isempty(Vnum)
            Vnum = Vnum(1);
            VsamSeq = Seq(1:VNJ(1));
            Vdel = VDJdata{j,Map.lDel(1)};
            
            %Choose the kappa/lampda DB
            if ~isempty(regexpi(Vname,'IGK')) 
                VrefSeq = DB.Vkmap{Vnum,1};
            else
                VrefSeq = DB.Vlmap{Vnum,1};
            end

            %Determine if there are extra nts left for padding
            ExtraLeft = VNJ(1)+ Vdel - length(VrefSeq);
            if ExtraLeft > 0
                VrefSeq = sprintf('%s%s',repmat('N',ExtraLeft,1),VrefSeq);
            elseif ExtraLeft < 0 %Delete some
                VrefSeq(1:abs(ExtraLeft)) = [];
            end

            %Get the deleted nts
            if Vdel > 0
                Vref3 = VrefSeq(end-Vdel+1:end);
                VrefSeq(end-Vdel+1:end) = [];
            else
                Vref3 = '';
            end

            %Assembe the V alignment string and save
            Valign = VrefSeq;
            Valign(VsamSeq == VrefSeq | VsamSeq == 'N') = '.';
            LightAlign{j,1} = sprintf('%s|%s',Valign,lower(Vref3));
        end

        %----------------------------------------------------------------------
        %Setting J alignment
        Jnum = VDJdata{j,Map.lGeneNum(2)};
        Jname = VDJdata{j,Map.lGeneName(2)};
        if ~isempty(Jnum)
            Jnum = Jnum(1);
            JsamSeq = Seq(end-VNJ(3)+1:end);
            JntDel = VDJdata{j,Map.lDel(2)};
            
            %Choose the kappa/lampda DB
            if ~isempty(regexpi(Jname,'IGK')) 
                JrefSeq = DB.Jkmap{Jnum,1};
            else
                JrefSeq = DB.Jlmap{Jnum,1};
            end
            
            %Determine if there are extra nts left for padding
            ExtraRight = VNJ(3) + JntDel - length(JrefSeq);
            if ExtraRight > 0
                JrefSeq = sprintf('%s%s',JrefSeq,repmat('N',ExtraRight,1));
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
            Jalign(JsamSeq == JrefSeq | JsamSeq == 'N') = '.';
            LightAlign{j,2} = sprintf('%s|%s',lower(Jref5),Jalign);
        end
    end
    
    LightAlignLoc = [Map.lValign Map.lJalign];
    VDJdata(:,LightAlignLoc) = LightAlign;
end

