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
%      'Update'        * 'y'       Update data such as CDR3, SHM, RefSeq
%                        'n'       Do not update other fields besides genes
%      'DJreserve'     *  9        Preserve 3' nts for D+J before V match
%      'Dreserve'      *  3        Perserve 5' nts for D before J match
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

function [VDJdata, BadLoc] = findVDJmatch(VDJdata, Map, DB, varargin)
BadLoc = repelem(false, size(VDJdata,1), 1);
if isempty(VDJdata); return; end

P = inputParser;
addParameter(P, 'Update', 'y', @(x) ismember(lower(x(1)), {'y', 'n'}));
addParameter(P, 'DJreserve', 9, @isnumeric);
addParameter(P, 'Dreserve', 3, @isnumeric);
parse(P, varargin{:});
Update = P.Results.Update;
DJreserve = P.Results.DJreserve; %Preserve last 9 for J.
Dreserve = P.Results.Dreserve;   %Preserve first 3 nts after V for D.

BadLoc = zeros(size(VDJdata, 1), 1, 'logical');
if ~contains(Map.Chain, 'H'); return; end

%Place VDJheader into separate variables due to parfor broadcast issues
SeqIdx      = Map.hSeq;
OverSeq5Idx = Map.hOverSeq5;
OverSeq3Idx = Map.hOverSeq3;
CDR3sIdx    = Map.hCDR3(3);
CDR3eIdx    = Map.hCDR3(4);
GeneNumIdx  = Map.hGeneNum;
GeneNameIdx = Map.hGeneName;
LengthIdx   = Map.hLength;
DelIdx      = Map.hDel;

%Place DB data into separate variables due to parfor broadcast issues.
M = getMapHeaderVar(DB.MapHeader);
AnchorIdx = M.Anchor;
Vmap = DB.Vmap;
Dmap = DB.Dmap;
Jmap = DB.Jmap;

%Begin finding the VDJ genes  
parfor j = 1:size(VDJdata, 1)
    MissRate = 0.15; %Start with 15% allowed mismatch rate. Place inside parfor to reset! 

    Tdata = VDJdata(j, :);
    Seq = Tdata{1, SeqIdx};
    if length(Seq) <= DJreserve %Sequence is too short
        BadLoc(j) = 1;
        continue
    end
    
    CDR3s = Tdata{1, CDR3sIdx};
    CDR3e = Tdata{1, CDR3eIdx};
    if isempty(CDR3s); CDR3s = 0; end
    if isempty(CDR3e); CDR3e = 0; end

    %Look for V gene for each seq
    Vnt = Seq(1:end-DJreserve); %Preserve nt to prevent matching V without D and J.
    CDR3s((CDR3s > length(Vnt)) | (CDR3s < 1)) = []; %Remove nonsensical locations
    Vmatch = findGeneMatch(Vnt, Vmap, 'V', MissRate, CDR3s);
    if isempty(Vmatch{1}) || max(Vmatch{1, 5}) <= 0 
        BadLoc(j) = 1;
        continue
    end
    Vlen = sum(Vmatch{1, 4}(1:2)); %Length of V segment

    %------------------------------------------------------------------
    %For short seq cutoff at 118W, if there's a possible C to W CDR3
    %region, then use the W position as a strong anchor for J alignment

    %Recalculate the CDR3s, 104C codon 1st nt
    V3del = Vmatch{3}(3); %V3' deletion
    VmapNum = Vmatch{1}(1);
    VrefCDR3 = Vmap{VmapNum, AnchorIdx}; %Location of C from right end of V
    CDR3s = Vlen + V3del - VrefCDR3 + 1;
    if CDR3s <= 0
        BadLoc(j) = 1;
        continue
    end

    %See if there's an in-frame 118W codon seed (if yes, use 'forceanchor')
    if max(CDR3e) == 0 %Find some CDR3e options
        CDR3e = strfind(Seq(Vlen+1:end), 'TGG') + 2 + Vlen;
        if isempty(CDR3e); CDR3e = 0; end
    end
    
    %CODING_NOTE: future release, try to remove this "forceanchor" as it will cause issues when dealing with out-of-frame junctions. 
    %Either skip it, or check BOTH force anchor + unforced anchor to determine which is the winning alignment.
    %2018-09-11
    
    ForceAnchor = '';
    if max(CDR3e) > 0
        %Determine in-frame TGGs
        InframeLoc = (mod(CDR3e-CDR3s-2, 3) == 0) & CDR3e > CDR3s;
        if any(InframeLoc)
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
        BadLoc(j) = 1;
        continue
    end
    Jlen = sum(Jmatch{4}(2:3));

    %Look for D gene for each seq
    Dnt = Seq(Vlen+1:end-Jlen);
    Dmatch = findGeneMatch(Dnt, Dmap, 'D', MissRate);
    if isempty(Dmatch{1})
        BadLoc(j) = 1;
        continue
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
        Over5Seq = [Tdata{1, OverSeq5Idx} Seq(1:VoverLen)];
        TrimSeq = 'y';
    else
        VoverLen = 0;
        Over5Seq = '';
    end
    JoverLen = JsamLMR(3) - JrefLMR(3);
    if JoverLen > 0
        JsamLMR(3) = JsamLMR(3) - JoverLen;
        Over3Seq = [Seq(end-JoverLen+1:end) Tdata{1, OverSeq3Idx}];
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
    Tdata(1, [SeqIdx OverSeq5Idx OverSeq3Idx]) = {Seq Over5Seq Over3Seq};
    Tdata(1, LengthIdx)   = num2cell(VMDNJ);
    Tdata(1, DelIdx)      = num2cell([V3del D5del D3del J5del]);
    Tdata(1, GeneNumIdx)  = [Vmatch(1, 1) Dmatch(1, 1) Jmatch(1, 1)];
    Tdata(1, GeneNameIdx) = [Vmatch(1, 2) Dmatch(1, 2) Jmatch(1, 2)];

    VDJdata(j, :) = Tdata;
end

if any(BadLoc)
    fprintf('%s: Bad sequences found.\n', mfilename);
    fprintf('  Removing Seq %d\n', find(BadLoc));
end

if strcmpi(Update(1), 'Y')
    KeepLoc = ~BadLoc;
    VDJdata(KeepLoc, :) = buildRefSeq(VDJdata(KeepLoc, :), Map, DB, 'H', 'germline', 'single');
    VDJdata(KeepLoc, :) = updateVDJdata(VDJdata(KeepLoc, :), Map, DB);
end