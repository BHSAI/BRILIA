%findVDJmatch will look for the best V, D, J gene match for a sequence in the
%following order: find V preserving nts for DJ, find J perserving nts for
%D, find D. This is good for an initial heavy chain annotation.
%
%  VDJdata = findVDJmatch(VDJdata, VDJheader, DB)
%
%  VDJdata = findVDJmatch(VDJdata, VDJheader, DB, Param, Value, ...)
%
%  [VDJdata, BadIdx] = findVDJmatch(...)
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
%      'DJreserve'     9           How may nts on right of Seq to preserve
%                                    for D+J gene matching
%      'Dreserve'      3           How may nts on right of Seq to preserve
%                                    for D gene matching
%
%  NOTE
%    If CDR3start and CDR3end anchor locations are provided in VDJdata, 
%    findVDJmatch will uses these location(s) to speed up V and J
%    alignments. Use seedCDR3position.m prior to this to benefit from
%    faster alignment. 
%
%    Using seeded alignment is faster BUT could cause this to miss the best
%    alignment. If you have the time, consider setting CDR3start and
%    CDR3end values to 0 or empty in VDJdata.
% 
%    After finding V, it will check to see if a 'TGG' of the 118W exists
%    after the 'TGT' of the conserved 104C, that is in frame. If it does
%    find one, it will force a J alignment such that the 118W are aligned.
%    This helps if the input seq is too short and ends exactly on the
%    conserved 118W. Otherwise, it either uses the seed from
%    seedCDR3position or a full J alignment to determine the J segment.
%
%  See also findVJmatch

function [VDJdata, varargout] = findVDJmatch(VDJdata, VDJheader, DB, varargin)
%Parse the input
P = inputParser;
addParameter(P, 'Update', 'y', @(x) ismember(lower(x(1)), {'y', 'n'}));
addParameter(P, 'DJreserve', 9, @isnumeric);
addParameter(P, 'Dreserve', 3, @isnumeric);
parse(P, varargin{:});
Update = P.Results.Update;
DJreserve = P.Results.DJreserve; %Preserve last 9 for J.
Dreserve = P.Results.Dreserve;   %Preserve first 3 nts after V for D.

BadIdx = zeros(size(VDJdata, 1), 1, 'logical');

%Place VDJheader into separate variables due to parfor broadcast issues
H = getHeavyHeaderVar(VDJheader);
if H.SeqLoc == 0 %Make sure header exists for right chain
    if nargout >=2
        varargout{1} = BadIdx;
    end
    return; 
end
SeqLoc      = H.SeqLoc;
OverSeq5Loc = H.OverSeq5Loc;
OverSeq3Loc = H.OverSeq3Loc;
CDR3sLoc    = H.CDR3Loc(3);
CDR3eLoc    = H.CDR3Loc(4);
GeneNumLoc  = H.GeneNumLoc;
GeneNameLoc = H.GeneNameLoc;
LengthLoc   = H.LengthLoc;
DelLoc      = H.DelLoc;

%Place DB data into separate variables due to parfor broadcast issues.
M = getMapHeaderVar(DB.MapHeader);
AnchorLoc = M.AnchorLoc;
Vmap = DB.Vmap;
Dmap = DB.Dmap;
Jmap = DB.Jmap;

%Begin finding the VDJ genes   
parfor j = 1:size(VDJdata, 1)
    %Start with 15% allowed mismatch rate. Place inside parfor to reset!
    MissRate = 0.15; 

    %Extract info from sliced VDJdata variable
    Tdata = VDJdata(j, :);
    Seq = Tdata{1, SeqLoc};
    if length(Seq) < DJreserve+1 %Sequence is too short
        BadIdx(j) = 1;
        continue; 
    end
    
    %Extract CDR3 start and end index. Set to 0 if empty.
    CDR3s = Tdata{1, CDR3sLoc};
    CDR3e = Tdata{1, CDR3eLoc};
    if isempty(CDR3s); CDR3s = 0; end
    if isempty(CDR3e); CDR3e = 0; end

    %Look for V gene for each seq
    Vnt = Seq(1:end-DJreserve); %Preserve nt to prevent matching V without D and J.
    CDR3s((CDR3s > length(Vnt)) | (CDR3s < 1)) = []; %Remove nonsensical locations
    Vmatch = findGeneMatch(Vnt, Vmap, 'V', MissRate, CDR3s);
    if isempty(Vmatch{1}) || max(Vmatch{1, 5}) <= 0 
        BadIdx(j) = 1;
        continue;
    end
    Vlen = sum(Vmatch{1, 4}(1:2)); %Length of V segment

    %------------------------------------------------------------------
    %For short seq cutoff at 118W, if there's a possible C to W CDR3
    %region, then use the W position as a strong anchor for J alignment

    %Recalculate the CDR3s, 104C codon 1st nt
    V3del = Vmatch{3}(3); %V3' deletion
    VmapNum = Vmatch{1}(1);
    VrefCDR3 = Vmap{VmapNum, AnchorLoc}; %Location of C from right end of V
    CDR3s = Vlen + V3del - VrefCDR3 + 1;
    if CDR3s <= 0 %Have an issue
        BadIdx(j) = 1;
        continue;
    end

    %See if there's an in-frame 118W codon seed (if yes, use 'forceanchor')
    if max(CDR3e) == 0 %Find some CDR3e options
        CDR3e = regexp(Seq(Vlen+1:end), 'TGG', 'end') + Vlen;
        if isempty(CDR3e); CDR3e = 0; end
    end
    ForceAnchor = '';
    if max(CDR3e) > 0
        %Determine in-frame TGGs
        InframeLoc = (mod(CDR3e-CDR3s-2, 3) == 0) & CDR3e > CDR3s;
        if max(InframeLoc) > 0
            CDR3e = CDR3e(InframeLoc);
            ForceAnchor = 'forceanchor'; %Ensure one of the anchor matches, UNLESS you get 0 matches.
        end
    end
    %End of ForceAnchor setting----------------------------------------

    MissRate = (Vlen - Vmatch{1, 4}(1) - Vmatch{1, 5}(1)) / (Vlen - Vmatch{1, 4}(1)); %Re-establish MissRate for D and J

    %Look for J gene for each seq
    Jnt = Seq(Vlen+Dreserve+1:end); %Preserve D NTs for later matching, and pad seq for OverhangMatch success.
    CDR3endTemp = CDR3e - Vlen - Dreserve; %Remember to shift CDR3e values
    CDR3endTemp((CDR3endTemp < 1) | (CDR3endTemp > length(Jnt))) = []; %Remove nonsensical locations
    Jmatch = findGeneMatch(Jnt, Jmap, 'J', MissRate, CDR3endTemp, ForceAnchor);
    if isempty(Jmatch{1}) || max(Jmatch{1, 5}) <= 0
        BadIdx(j) = 1;
        continue;
    end
    Jlen = sum(Jmatch{4}(2:3));

    %Look for D gene for each seq
    Dnt = Seq(Vlen+1:end-Jlen);
    Dmatch = findGeneMatch(Dnt, Dmap, 'D', MissRate);
    if isempty(Dmatch{1})
        BadIdx(j) = 1;
        continue;
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
        Over5Seq = [Tdata{1, OverSeq5Loc} Seq(1:VoverLen)];
        TrimSeq = 'y';
    else
        VoverLen = 0;
        Over5Seq = '';
    end
    JoverLen = JsamLMR(3) - JrefLMR(3);
    if JoverLen > 0
        JsamLMR(3) = JsamLMR(3) - JoverLen;
        Over3Seq = [Seq(end-JoverLen+1:end) Tdata{1, OverSeq3Loc}];
        TrimSeq = 'y';
    else
        JoverLen = 0;
        Over3Seq = '';
    end
    if TrimSeq == 'y'
        Seq = Seq(VoverLen+1:end-JoverLen);
    end
    
    %Extract VMDNJ lengths
    Vlen = sum(VsamLMR(1:2));
    Jlen = sum(JsamLMR(2:3));
    Dlmr = cell2mat(Dmatch(1, 4)); %J left mid right segment
    VMDNJ = [Vlen Dlmr Jlen]; %V, Nvd, D, Ndj, J lengths

    %Extract germline deletion info
    V3del = Vmatch{1, 3}(3);
    D5del = Dmatch{1, 3}(1);
    D3del = Dmatch{1, 3}(3);
    J5del = Jmatch{1, 3}(1);

    %Extract the gene family map number and family resolution
    Tdata(1, [SeqLoc OverSeq5Loc OverSeq3Loc]) = {Seq Over5Seq Over3Seq};
    Tdata(1, LengthLoc)   = num2cell(VMDNJ);
    Tdata(1, DelLoc)      = num2cell([V3del D5del D3del J5del]);
    Tdata(1, GeneNumLoc)  = [Vmatch(1, 1) Dmatch(1, 1) Jmatch(1, 1)];
    Tdata(1, GeneNameLoc) = [Vmatch(1, 2) Dmatch(1, 2) Jmatch(1, 2)];

    VDJdata(j, :) = Tdata;
end

%If there are errors, show them now.
BadLoc = find(BadIdx == 1);
for b = 1:length(BadLoc)
    fprintf('%s: Bad sequence # %d\n', mfilename, BadLoc(b));
end

%Update VDJdata 
if upper(Update(1)) == 'Y'
    UpdateIdx = ~BadIdx;
    VDJdata(UpdateIdx, :) = buildRefSeq(VDJdata(UpdateIdx, :), VDJheader, DB, 'H', 'germline', 'single');
    VDJdata(UpdateIdx, :) = updateVDJdata(VDJdata(UpdateIdx, :), VDJheader, DB);
end

if nargout >=2
    varargout{1} = BadIdx;
end
