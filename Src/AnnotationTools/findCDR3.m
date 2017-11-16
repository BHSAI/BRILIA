%findCDR3 will return the CDR3 information, such as AA, AA length, and
%start and end location of this region for the Sequence in VDJdata.
%
%  VDJdata = findCDR3(VDJdata, VDJheader, DB)
%
%  VDJdata = findCDR3(VDJdata, VDJheader, DB, Option)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    VDJheader: main BRILIA header cell
%    DB: Gene database structure (getCurrentDatabase.m)
%    Option ['IMGT','anchor]: specify whether to return the IMGT's
%      definition of CDR3 or the CDR3 anchor. The anchor positions include
%      the 104C and 118F/W codons. Default is 'IMGT'.
%
%  OUTPUT
%    VDJdata: modified VDJdata with the CDR3 information filled
%
%  NOTE
%    This must be run after finding V(D)J annotations.
%
%  See also findCDR1, findCDR2

function VDJdata = findCDR3(VDJdata,Map,DB,varargin)
%Determine which CDR3 output to return - IMGT's CDR3 or CDR3 anchor
Option = 'imgt';
if ~isempty(varargin) && ischar(varargin{1})
    Option = lower(varargin{1});
end

%Determine chain and extract key locations
M = getMapHeaderVar(DB.MapHeader);

for k = 1:length(Map.Chain)
    %Decide what header locator to use
    if Map.Chain(k) == 'H'
        DelLoc  = Map.hDel;
        GeneNumLoc  = Map.hGeneNum; 
        SeqLoc  = Map.hSeq;
        LengthLoc  = Map.hLength;
        CDR3Loc = Map.hCDR3;
    else
        DelLoc  = Map.lDel;
        GeneNumLoc  = Map.lGeneNum; 
        SeqLoc  = Map.lSeq;
        LengthLoc  = Map.lLength;
        CDR3Loc = Map.lCDR3;
    end
    
    for j = 1:size(VDJdata,1)
        %Get basic annotation information
        Seq = VDJdata{j,SeqLoc};
        Vlen = VDJdata{j,LengthLoc(1)};
        Jlen = VDJdata{j,LengthLoc(end)};
        Vdel = VDJdata{j,DelLoc(1)};
        Jdel = VDJdata{j,DelLoc(end)};
        VmapNum = VDJdata{j,GeneNumLoc(1)};
        JmapNum = VDJdata{j,GeneNumLoc(end)};

        %Make sure all info is provided for CDR3 calculation 
        if isempty(Seq); continue; end
        if isempty(Vlen); continue; end
        if isempty(Jlen); continue; end
        if isempty(Vdel); continue; end
        if isempty(Jdel); continue; end
        if isempty(VmapNum); continue; end
        if isempty(JmapNum); continue; end
        
        if Map.Chain(k) == 'H' %Heavy chain, use Vmap Jmap
            VrefCDR3 = DB.Vmap{VmapNum(1),M.AnchorLoc}; %Location of C from right end of V
            JrefCDR3 = DB.Jmap{JmapNum(1),M.AnchorLoc}; %Location of W from left end of J
        else %Light chain, determine locus and u se Vkmap, Vlmap, Jkmap, Jlmap
            %Determine the locus and anchor loc
            Vname = VDJdata{j,Map.lGeneName(1)};
            if ~isempty(regexpi(Vname,'IGKV','once')) %Kappa
                VrefCDR3 = DB.Vkmap{VmapNum(1),M.AnchorLoc}; %Location of C from right end of V
                JrefCDR3 = DB.Jkmap{JmapNum(1),M.AnchorLoc}; %Location of F from left end of J
            else %Lambda
                VrefCDR3 = DB.Vlmap{VmapNum(1),M.AnchorLoc}; %Location of C from right end of V
                JrefCDR3 = DB.Jlmap{JmapNum(1),M.AnchorLoc}; %Location of F from left end of J
            end
        end
        
        %Recalculate the CDR3start, 104C codon 1st nt
        CDR3s = Vlen + Vdel - VrefCDR3 + 1;
        if CDR3s < 0; CDR3s = 0; end

        %Recalculate the CDR3end, 118W/F codon 3rd nt
        CDR3e = length(Seq) - (Jlen + Jdel) + JrefCDR3 + 2;
        if CDR3e < 0; CDR3e = 0; end

        %Make sure CDR3 has >= 5 aa or 15 nt, and does not have 0 value.
        CDR3aa = '';
        if CDR3e > 0 && CDR3s > 0 && (CDR3e-CDR3s+1) > 15
            if CDR3e > length(Seq) %Just in case you have a cutoff seq
                CDR3nt = Seq(CDR3s:end);
            else
                CDR3nt = Seq(CDR3s:CDR3e);
            end
            CDR3aa = convNT2AA(CDR3nt,'ACGTonly',false);
        end
        
        %If the default option of IMGT is selected, need to shift positions
        if strcmpi(Option,'imgt')
            CDR3aa = CDR3aa(2:end-1);
            CDR3s = CDR3s+3;
            CDR3e = CDR3e-3;
        end

        %Fill in the VDJdata
        VDJdata(j,CDR3Loc) = {CDR3aa length(CDR3aa) CDR3s CDR3e};
    end
end
