%buildVDJalignment will return the VDJ alignment results for each sequence.
%An "-" is a match with the germline gene, a lower case are deleted nts
%from the germline gene, and a capital letter is a mismatch with the
%germline gene, and an "|" marks where the germline gene was deleted.
%
%  VDJdata = buildVDJalignment(VDJdata, NewHeader)
%
%  VDJdata = buildVDJalignment(VDJdata, NewHeader, Vmap, Dmap, Jmap)

function VDJdata = buildVDJalignment(VDJdata,NewHeader,varargin)
%Extract the VDJ database
if length(varargin) >= 3
    [Vmap,Dmap,Jmap] = deal(varargin{1:3});
else
    [Vmap,Dmap,Jmap] = getCurrentDatabase;
end
getHeaderVar;

for j = 1:size(VDJdata,1)
    try
        Seq = VDJdata{j,SeqLoc};
        VMDNJ = cell2mat(VDJdata(j,LengthLoc));

        %----------------------------------------------------------------------
        %Setting V alignment
        Vnum = VDJdata{j,FamNumLoc(1)}(1);
        if ~isempty(Vnum)
            VsamSeq = Seq(1:VMDNJ(1));
            VrefSeq = Vmap{Vnum,1};
            VntDel = VDJdata{j,DelLoc(1)};

            %Determine if there are extra nts left for padding
            ExtraLeft = VMDNJ(1)+ VntDel - length(VrefSeq);
            if ExtraLeft > 0
                VrefSeq = sprintf('%s%s',repmat('X',ExtraLeft,1),VrefSeq);
            elseif ExtraLeft < 0 %Delete some
                VrefSeq(1:abs(ExtraLeft)) = [];
            end

            if VntDel > 0
                Vref3 = VrefSeq(end-VntDel+1:end);
                VrefSeq(end-VntDel+1:end) = [];
            else
                Vref3 = '';
            end

            Valign = VrefSeq;
            Valign(VsamSeq == VrefSeq | VsamSeq == 'X') = '-';
            Valign = [Valign '|' lower(Vref3)];
            VDJdata{j,AlignLoc(1)} = Valign;
        end

        %----------------------------------------------------------------------
        %Setting D alignment
        Dnum = VDJdata{j,FamNumLoc(2)}(1);
        if ~isempty(Dnum)
            DsamSeq = Seq(sum(VMDNJ(1:2))+1:sum(VMDNJ(1:3)));

            %Extract the germline seq
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
            Dalign(DsamSeq==DrefSeq | DsamSeq == 'X') = '-';
            Dalign = [lower(Dref5) '|' Dalign '|' lower(Dref3)];
            VDJdata{j,AlignLoc(2)} = Dalign;
        end

        %----------------------------------------------------------------------
        %Setting J alignment
        Jnum = VDJdata{j,FamNumLoc(3)}(1);
        if ~isempty(Jnum)
            JsamSeq = Seq(end-VMDNJ(5)+1:end);
            JrefSeq = Jmap{Jnum,1};
            JntDel = VDJdata{j,DelLoc(4)};

            %Determine if there are extra nts left for padding
            ExtraRight = VMDNJ(5) + JntDel - length(JrefSeq);
            if ExtraRight > 0
                JrefSeq = sprintf('%s%s',JrefSeq,repmat('X',ExtraRight,1));
            elseif ExtraRight < 0
                JrefSeq(end-abs(ExtraRight)+1:end)= [];
            end

            if JntDel > 0
                Jref5 = JrefSeq(1:JntDel);
                JrefSeq(1:JntDel) = [];
            else
                Jref5 = '';
            end

            Jalign = JrefSeq;
            Jalign(JsamSeq == JrefSeq | JsamSeq == 'X') = '-';
            Jalign = [lower(Jref5) '|' Jalign];
            VDJdata{j,AlignLoc(3)} = Jalign;
        end
    catch
        ErrorMsg = sprintf('Errored at %s, sequence # %d',mfilename,j);
        disp(ErrorMsg);
        VDJdata{j,MiscLoc} = ErrorMsg;
        pause
    end
end
