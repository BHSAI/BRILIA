%findVJmatch will look for the best V and J gene match for a sequence in
%the following order: find V preserving nts for J, and then find J This is
%good for an initial light chain annotation.
%
%  VDJdata = findVJmatch(VDJdata,VDJheader,DB)
%
%  VDJdata = findVJmatch(VDJdata,VDJheader,DB,Param,Value,...)
%
%  [VDJdata,BadIdx] = findVJmatch(...)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    VDJheader: main BRILIA header cell
%    DB: Gene database structure (getCurrentDatabase.m)
%
%      Param           Value       Description
%      --------------- ---------   ----------------------------------------
%      'Update'        'y' 'n'     Yes or No to updating other fields
%                                    besides the gene fields, such as CDR3,
%                                    SHM, RefSeq.
%      'Jreserve'      9           How may nts on right of Seq to preserve
%                                    for J gene matching
%
%  NOTE
%    If CDR3start and CDR3end anchor locations are provided in VDJdata,
%    findVJmatch will uses these location(s) to speed up V and J
%    alignments. Use seedCDR3position.m prior to this to benefit from
%    faster alignment. 
%
%    Using seeded alignment is faster BUT could cause this to miss the best
%    alignment. If you have the time, consider setting CDR3start and
%    CDR3end values to 0 or empty in VDJdata.
% 
%    After finding V, it will check to see if a 'TTT' or 'TTC' of the 118F
%    exists after the 'TGT' of the conserved 104C, that is in frame. If it
%    does find one, it will force a J alignment such that the 118F are
%    aligned. This helps if the input seq is too short and ends exactly on
%    the conserved 118F. Otherwise, it either uses the seed from
%    seedCDR3position or a full J alignment to determine the J segment.
%
%  See also findVDJmatch

function [VDJdata,varargout] = findVJmatch(VDJdata,Map,DB,varargin)
%Parse the input
P = inputParser;
addParameter(P,'Update','y',@(x) ismember(lower(x(1)),{'y','n'}));
addParameter(P,'Jreserve',9,@isnumeric);
parse(P,varargin{:});
Update = P.Results.Update;
Jreserve = P.Results.Jreserve;

BadIdx = zeros(size(VDJdata,1),1,'logical');

%Find headers since parfor can't handle it
if Map.lSeq == 0 %Make sure header exists for right chain   
    if nargout >=2
        varargout{1} = BadIdx;
    end
    return; 
end
SeqLoc      = Map.lSeq;
OverSeq5Loc = Map.lOverSeq5;
OverSeq3Loc = Map.lOverSeq3;
CDR3sLoc    = Map.lCDR3(3);
CDR3eLoc    = Map.lCDR3(4);
GeneNumLoc  = Map.lGeneNum;
GeneNameLoc = Map.lGeneName;
LengthLoc   = Map.lLength;
DelLoc      = Map.lDel;

%Place DB into separate tables due to broadcast issues.
%Find headers since parfor can't handle it
M = getMapHeaderVar(DB.MapHeader);
AnchorLoc   = M.Anchor;
Vxmap = [DB.Vkmap; DB.Vlmap]; %Combined for simplicity
VkCount = size(DB.Vkmap,1);   %Marks how many are for Vk vs. Vl
Jkmap = DB.Jkmap;
Jlmap = DB.Jlmap;

%Begin finding the VJ genes   
parfor j = 1:size(VDJdata,1)
    %Extract info from sliced VDJdata variable
    Tdata = VDJdata(j,:);
    Seq = Tdata{1,SeqLoc};
    if length(Seq) < Jreserve+1 %Seq is too short
        BadIdx(j) = 1;
        continue; 
    end
    CDR3s = Tdata{1,CDR3sLoc};
    CDR3e = Tdata{1,CDR3eLoc};
    if isempty(CDR3s); CDR3s = 0; end
    if isempty(CDR3e); CDR3e = 0; end

    MissRate = 0.15; %Start with 15% miss rate for the V segment. Place inside parfor loop to reset!

    %Look for V gene for each seq
    Vnt = Seq(1:end-Jreserve); %Preserve nt to prevent matching V without J.
    CDR3s((CDR3s > length(Vnt)) | (CDR3s < 1)) = []; %Remove nonsensical locations
    Vmatch = findGeneMatch(Vnt,Vxmap,'V',MissRate,CDR3s);
    Vlen = sum(Vmatch{1,4}(1:2)); %Length of V segment
    
    %------------------------------------------------------------------
    %For short seq cutoff at 118F, if there's a possible C to F CDR3
    %region, then use the F position as a strong anchor for J alignment

    %Recalculate the CDR3start, 104C codon 1st nt
    Vdel = Vmatch{3}(3); %V3' deletion
    VmapNum = Vmatch{1}(1);
    VrefCDR3 = Vxmap{VmapNum,AnchorLoc}; %Location of C from right end of V
    CDR3s = Vlen + Vdel - VrefCDR3 + 1;
    if CDR3s <= 0 %Have an issue
        BadIdx(j) = 1;
        continue;
    end

    %See if there's an in-frame 118F codon seed (if yes, use 'forceanchor')
    if max(CDR3e) == 0 %Find some CDR3end options
        CDR3e = regexp(Seq(Vlen+1:end),'TT[TC]','end') + Vlen;
        if isempty(CDR3e); CDR3e = 0; end
    end
    ForceAnchor = '';
    if max(CDR3e) > 0
        %Determine in-frame TGGs
        InframeLoc = (mod(CDR3e-CDR3s-2,3) == 0) & CDR3e > CDR3s;
        if max(InframeLoc) > 0
            CDR3e = CDR3e(InframeLoc);
            ForceAnchor = 'forceanchor'; %Ensure once of the anchor matches, UNLESS you get 0 matches.
        end
    end
    %End of ForceAnchor setting----------------------------------------
    
    %Determine actual miss rate and also locus
    MissRate = (Vlen - Vmatch{1,4}(1) - Vmatch{1,5}(1)) / (Vlen - Vmatch{1,4}(1));
    GeneName = Vmatch{1,2};
    if ~isempty(regexpi(GeneName,'IGKV')) %Kappa
        Locus = 'Kappa';
    else %Lambda
        Locus = 'Lambda';
        Vmatch{1,1} = Vmatch{1,1} - VkCount; %Renumber relative to Vlmap
    end
    
    %Look for J gene for each seq
    Jnt = Seq(Vlen+1:end); %Preserve D NTs for later matching, and pad seq for OverhangMatch success.
    CDR3endTemp = CDR3e - Vlen; %Remember to shift CDR3end values
    CDR3endTemp((CDR3endTemp < 1) | (CDR3endTemp > length(Jnt))) = []; %Remove nonsensical locations
    if Locus(1) == 'K'
        Jmatch = findGeneMatch(Jnt,Jkmap,'J',MissRate,CDR3endTemp,ForceAnchor);
    else
        Jmatch = findGeneMatch(Jnt,Jlmap,'J',MissRate,CDR3endTemp,ForceAnchor);
    end
    
    %Begin trimming sequences that go beyond the V or J genes
    VsamLMR = Vmatch{4};
    VrefLMR = Vmatch{3};
    JsamLMR = Jmatch{4};
    JrefLMR = Jmatch{3};
    VoverLen = VsamLMR(1) - VrefLMR(1);
    TrimSeq = 'n';
    if VoverLen > 0
        VsamLMR(1) = VsamLMR(1) - VoverLen;
        Over5Seq = [Tdata{1,OverSeq5Loc} Seq(1:VoverLen)];
        TrimSeq = 'y';
    else
        VoverLen = 0;
        Over5Seq = '';
    end
    JoverLen = JsamLMR(3) - JrefLMR(3);
    if JoverLen > 0
        JsamLMR(3) = JsamLMR(3) - JoverLen;
        Over3Seq = [Seq(end-JoverLen+1:end) Tdata{1,OverSeq3Loc}];
        TrimSeq = 'y';
    else
        JoverLen = 0;
        Over3Seq = '';
    end
    if TrimSeq == 'y'
        Seq = Seq(VoverLen+1:end-JoverLen);
    end
    
    %Extract VNJ lengths, excluding any overhanging sequences
    Vlen = sum(VsamLMR(1:2));
    Jlen = sum(JsamLMR(2:3));
    Nlen = length(Seq) - Vlen - Jlen;
    VNJ = [Vlen Nlen Jlen]; %V, N, J lengths

    %Extract germline deletion info
    Vdel = Vmatch{1,3}(3);
    Jdel = Jmatch{1,3}(1);

    %Extract update the fields in Tdata and save to VDJdata
    Tdata(1,[SeqLoc OverSeq5Loc OverSeq3Loc]) = {Seq Over5Seq Over3Seq};
    Tdata(1,LengthLoc)   = num2cell(VNJ);
    Tdata(1,DelLoc)      = num2cell([Vdel Jdel]);
    Tdata(1,GeneNumLoc)  = [Vmatch(1,1) Jmatch(1,1)];
    Tdata(1,GeneNameLoc) = [Vmatch(1,2) Jmatch(1,2)];

    VDJdata(j,:) = Tdata;
end

%If there are errors, show them now.
BadLoc = find(BadIdx == 1);
for b = 1:length(BadLoc)
    fprintf('%s: Bad sequence # %d\n',mfilename,BadLoc(b));
end

%Update VDJdata 
if upper(Update(1)) == 'Y'
    UpdateIdx = ~BadIdx;
    VDJdata(UpdateIdx,:) = buildRefSeq(VDJdata(UpdateIdx,:),Map,DB,'L','germline','single');
    VDJdata(UpdateIdx,:) = updateVDJdata(VDJdata(UpdateIdx,:),Map,DB);    
end

if nargout >=2
    varargout{1} = BadIdx;
end
