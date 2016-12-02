%findCDR3 will search for the CDR3 regions, and fill in the CDR3 info into
%VDJdata, such as CDR3 nt length, Start of CDR3 104C codon, End of CDR3 118
%codon.
%
%  VDJdata = findCDR3(VDJdata,NewHeader)
function VDJdata = findCDR3(VDJdata,NewHeader,varargin)
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
    
    %Determine the reference gene location of 104C
    if ischar(VDJdata{j,FamNumLoc(1)})
        VmapNum = eval(VDJdata{j,FamNumLoc(1)}); %V map number
    else
        VmapNum = VDJdata{j,FamNumLoc(1)}; %V map number
    end
    VmapNum = VmapNum(1);
    VrefCDR3 = Vmap{VmapNum(1),end}; %Location of C from right end of V
    Vdel = VDJdata{j,DelLoc(1)};
    Vlen = VDJdata{j,LengthLoc(1)};
    CDR3start = Vlen + Vdel - VrefCDR3 + 1;
    
    %Determine the reference gene location of 118W.
    if ischar(VDJdata{j,FamNumLoc(end)})
        JmapNum = eval(VDJdata{j,FamNumLoc(end)});
    else
        JmapNum = VDJdata{j,FamNumLoc(end)};
    end
    JmapNum = JmapNum(1);
    JrefCDR3 = Jmap{JmapNum(1),end}; %Location of F/W from left end of J
    Jdel = VDJdata{j,DelLoc(end)};
    Jlen = VDJdata{j,LengthLoc(end)};
    CDR3end = length(Seq) - (Jlen + Jdel) + JrefCDR3 + 2;
    
    if (CDR3start < 1) || CDR3end < sum(VMDNJ(1:4))
        VDJdata{j,CDR3Loc(1)} = '';
        VDJdata{j,CDR3Loc(2)} = 0;
        VDJdata{j,CDR3Loc(3)} = 0;
        VDJdata{j,CDR3Loc(4)} = 0;        
        continue %CDR3 is not defined well
    else
        if CDR3end > length(Seq)
            CDR3nt = Seq(CDR3start:end);
        else
            CDR3nt = Seq(CDR3start:CDR3end);
        end
        CDR3aa = nt2aa(CDR3nt,'ACGTonly','false');

        %Fill in the VDJdata
        VDJdata{j,CDR3Loc(1)} = CDR3aa;
        VDJdata{j,CDR3Loc(2)} = length(CDR3aa);
        VDJdata{j,CDR3Loc(3)} = CDR3start;
        VDJdata{j,CDR3Loc(4)} = CDR3end;        
    end
end