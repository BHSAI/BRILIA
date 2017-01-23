%findCDR3 will recheck the locations of the CDR3 region, translate the CDR3
%nt into aa if possible, and reset CDR3 start and end locations based on
%the 104C and 118W codons expected by the reference sequence annotations.
%This must be run after finding VDJ annotations.
%
%  VDJdata = findCDR3(VDJdata, NewHeader)
%
%  VDJdata = findCDR3(VDJdata, NewHeader, Vmap, Dmap, Jmap)
%
%  [VDJdata, BadIdx] = findCDR3(VDJdata, NewHeader)

function [VDJdata,varargout] = findCDR3(VDJdata,NewHeader,varargin)
%Extract the VDJ database
if length(varargin) >= 3
    [Vmap,Dmap,Jmap] = deal(varargin{1:3});
else
    [Vmap,Dmap,Jmap] = getCurrentDatabase;
end
getHeaderVar;

BadIdx = zeros(size(VDJdata,1),'logical');
for j = 1:size(VDJdata,1)
    try
        %Get basic annotation information
        Seq = VDJdata{j,SeqLoc};
        VMDNJ = cell2mat(VDJdata(j,LengthLoc));
        Vlen = VMDNJ(1);
        Vdel = VDJdata{j,DelLoc(1)};
        Jlen = VMDNJ(5);
        Jdel = VDJdata{j,DelLoc(4)};
        
        %Recalculate the CDR3start, 104C codon 1st nt
        VmapNum = VDJdata{j,FamNumLoc(1)}(1);
        VrefCDR3 = Vmap{VmapNum,10}; %Location of C from right end of V
        CDR3start = Vlen + Vdel - VrefCDR3 + 1;
        if CDR3start < 0; CDR3start = 0; end
        
        %Recalculate the CDR3end, 118W/F codon 3rd nt
        JmapNum = VDJdata{j,FamNumLoc(end)}(1);
        JrefCDR3 = Jmap{JmapNum,10}; %Location of F/W from left end of J
        CDR3end = length(Seq) - (Jlen + Jdel) + JrefCDR3 + 2;
        if CDR3end < 0; CDR3end = 0; end
        
        %Make sure CDR3 has >= 4 aa, and does not have 0 value.
        CDR3aa = '';
        if (CDR3end * CDR3start) > 0 && (CDR3end-CDR3start > 11)
            if CDR3end > length(Seq) %Just in case you have a cutoff seq
                CDR3nt = Seq(CDR3start:end);
            else
                CDR3nt = Seq(CDR3start:CDR3end);
            end

            %Replace any ambiguous characters as N
            AmbigLoc = regexpi(CDR3nt,'[^ACGTU]');
            CDR3nt(AmbigLoc) = 'N';
            CDR3aa = nt2aa(CDR3nt,'ACGTonly','false');
        end
        
        %Fill in the VDJdata
        VDJdata(j,CDR3Loc) = {CDR3aa length(CDR3aa) CDR3start CDR3end};
    catch
        ErrorMsg = sprintf('Errored at %s, sequence # %d',mfilename,j);
        disp(ErrorMsg);
        VDJdata{j,MiscLoc} = ErrorMsg;
        BadIdx(j) = 1;   
    end
end

if nargout >=2
    varargout{1} = BadIdx;
end