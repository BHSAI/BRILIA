%findVDJmatch will look for the best V,D,J gene match for a sequence in the
%following order: find V preserving nts for DJ, find J perserving nts for
%D, find D. This is good for getting an initial guess VDJ annotation.
%
%  VDJdata = findVDJmatch(VDJdata,VDJheader)
%
%  VDJdata = findVDJmatch(VDJdata,VDJheader,Vmap,Dmap,Jmap)
%
%  VDJdata = findVDJmatch(VDJdata,VDJheader,Vmap,Dmap,Jmap,'update')
%
%  [VDJdata,BadIdx] = findVDJmatch(...)
%
%  INPUT
%    'update': after finding VDJ matches, will assemble the reference
%      germline sequence and also fill in other fields of VDJdata such as
%      SHM, alignment, CDR3, etc. Without update, it will only fill in
%      fields related to VDJ genes.
%
%  NOTE
%    If CDR3start and CDR3end locations are provided in VDJdata,
%    findVDJmatch will uses these location(s) to speed up V and J
%    alignments. Use seedCDR3position.m prior to this to benefit from
%    faster alignment. 
%
%    Using seeded alignment is faster BUT could cause this to miss best
%    alignment. If you have the time, consider setting CDR3start and
%    CDR3end values to 0 or empty in VDJdata.
% 
%    After finding V, it will check to see if a 'TGG' of the 118W exists
%    after the 'TGT' of the conserved 104C, that is in frame. If it does
%    find one, it will force a J alignment such that the 118W are aligned.
%    This helps if the input seq is too short and ends exactly on the
%    conserved 118W. Otherwise, it either uses the seed from
%    seedCDR3positoin or a full J alignment to determine the J segment.
%
function [VDJdata,varargout] = findVDJmatch(VDJdata,VDJheader,varargin)
%Look if update was specified first
UpdateLoc = findCell(varargin,'update','MatchCase','any');
if UpdateLoc(1) ~= 0 %There is an update input
    Update = 'y';
    varargin(UpdateLoc) = [];
else
    Update = 'n';
end

%Parse the input
P = inputParser;
addOptional(P,'Vmap',{},@iscell);
addOptional(P,'Dmap',{},@iscell);
addOptional(P,'Jmap',{},@iscell);
addParameter(P,'DJreserve',9,@isnumeric);
addParameter(P,'Dreserve',3,@isnumeric);
parse(P,varargin{:});
Vmap = P.Results.Vmap; 
Dmap = P.Results.Dmap;
Jmap = P.Results.Jmap;
DJreserve = P.Results.DJreserve; %Preserve last 9 for J.
Dreserve = P.Results.Dreserve;   %Preserve first 3 nts after V for D.

%Extract the VDJ database, if none is specified.
if isempty(Vmap) || isempty(Dmap) || isempty(Jmap)
    [Vmap, Dmap, Jmap] = getCurrentDatabase('change');
    [Vmap, Dmap, Jmap] = filterRefGene(Vmap,Dmap,Jmap);
end

%Bring some H = getHeaderVar(VDJheader); variables out, since parfor can't handle it.
SeqLoc = findCell(VDJheader,{'nucleotide','Seq'});
LengthLoc = findCell(VDJheader,{'Length_V','Length_Nvd','Length_D','Length_Ndj','Length_J'});
DelLoc = findCell(VDJheader,{'V_Deletion3','D_Deletion5','D_Deletion3','J_Deletion5'});
FamNumLoc = findCell(VDJheader,{'V_MapNum','D_MapNum','J_MapNum'});
FamLoc = findCell(VDJheader,{'V_GeneName','D_GeneName','J_GeneName'});
MiscLoc = findCell(VDJheader,{'Misc'});
CDR3startLoc = findCell(VDJheader,{'CDR3_Start'});
CDR3endLoc = findCell(VDJheader,{'CDR3_End'});

%Begin finding the VDJ genes   
BadIdx = zeros(size(VDJdata,1),1,'logical');
parfor j = 1:size(VDJdata,1)
    %Extract info from sliced VDJdata variable
    Tdata = VDJdata(j,:);
    Seq = Tdata{1,SeqLoc};
    CDR3start = Tdata{1,CDR3startLoc};
    if isempty(CDR3start); CDR3start = 0; end
    CDR3end = Tdata{1,CDR3endLoc};
    if isempty(CDR3end); CDR3end = 0; end

    MissRate = 0.15; %Start with 15% miss rate for the V segment. Place inside parfor loop to reset!

    %Look for V gene for each seq
    Vnt = Seq(1:end-DJreserve); %Preserve nt to prevent matching V without D and J.
    AllowedMiss = ceil(MissRate * (length(Vnt))); %Number of pt mutations allowed
    CDR3start((CDR3start > length(Vnt)) | (CDR3start < 1)) = []; %Remove nonsensical locations
    Vmatch = findGeneMatch(Vnt,Vmap,'V',AllowedMiss,CDR3start);
    Vlen = sum(Vmatch{4}(1:2)); %Length of V segment

    %------------------------------------------------------------------
    %For short seq cutoff at 118W, if there's a possible C to W CDR3
    %region, then use the W position as a strong anchor for J alignment

    %Recalculate the CDR3start, 104C codon 1st nt
    Vdel = Vmatch{3}(3); %V3' deletion
    VmapNum = Vmatch{1}(1);
    VrefCDR3 = Vmap{VmapNum,10}; %Location of C from right end of V
    CDR3start = Vlen + Vdel - VrefCDR3 + 1;
    if CDR3start <= 0 %Have an issue
        BadIdx(j) = 1;
        continue;
    end

    %See if there's an in-frame 118W codon seed (if yes, use 'forceanchor')
    if max(CDR3end) == 0 %Find some CDR3end options
        CDR3end = regexp(Seq(Vlen+1:end),'TGG','end') + Vlen;
        if isempty(CDR3end); CDR3end = 0; end
    end
    ForceAnchor = '';
    if max(CDR3end) > 0
        %Determine in-frame TGGs
        InframeLoc = (mod(CDR3end-CDR3start-2,3) == 0) & CDR3end > CDR3start;
        if max(InframeLoc) > 0
            CDR3end = CDR3end(InframeLoc);
            ForceAnchor = 'forceanchor'; %Ensure once of the anchor matches, UNLESS you get 0 matches.
        end
    end
    %End of ForceAnchor setting----------------------------------------

    MissRate = (Vlen - Vmatch{1,5}(1))/Vlen; %Re-establish MissRate for D and J

    %Look for J gene for each seq
    Jnt = Seq(Vlen+Dreserve+1:end); %Preserve D NTs for later matching, and pad seq for OverhangMatch success.
    AllowedMiss = ceil(MissRate * length(Jnt)); %Number of pt mutations allowed
    CDR3endTemp = CDR3end - Vlen - Dreserve; %Remember to shift CDR3end values
    CDR3endTemp((CDR3endTemp < 1) | (CDR3endTemp > length(Jnt))) = []; %Remove nonsensical locations
    Jmatch = findGeneMatch(Jnt,Jmap,'J',AllowedMiss,CDR3endTemp,ForceAnchor);
    Jlen = sum(Jmatch{4}(2:3));

    %Look for D gene for each seq
    Dnt = Seq(Vlen+1:end-Jlen);
    AllowedMiss = ceil(MissRate * length(Dnt));
    Dmatch = findGeneMatch(Dnt,Dmap,'D',AllowedMiss);

    %Extract VMDNJ lengths
    Vlmr = cell2mat(Vmatch(1,4)); %V left mid right segment
    Jlmr = cell2mat(Jmatch(1,4)); %D left mid right segment
    Dlmr = cell2mat(Dmatch(1,4)); %J left mid right segment
    VMDNJ = [sum(Vlmr(1:2)) Dlmr sum(Jlmr(2:3))]; %V, Nvd, D, Ndj, J lengths
    Tdata(1,LengthLoc) = num2cell(VMDNJ);

    %Extract germline deletion info
    Vdel = Vmatch{1,3}(3);
    D5del = Dmatch{1,3}(1);
    D3del = Dmatch{1,3}(3);
    Jdel = Jmatch{1,3}(1);
    Tdata(1,DelLoc) = num2cell([Vdel D5del D3del Jdel]);

    %Extract the gene family map number and family resolution
    Tdata(1,FamNumLoc) = [Vmatch(1,1) Dmatch(1,1) Jmatch(1,1)];
    Tdata(1,FamLoc) = [Vmatch(1,2) Dmatch(1,2) Jmatch(1,2)];

    VDJdata(j,:) = Tdata;
end

%If there are errors, show them now.
BadLoc = find(BadIdx == 1);
for b = 1:length(BadLoc)
    ErrorMsg = sprintf('Errored at %s, sequence # %d',mfilename,BadLoc(b));
    disp(ErrorMsg);
    VDJdata{BadLoc(b),MiscLoc} = ErrorMsg;
end

%Update VDJdata 
if Update == 'y'
    UpdateIdx = ~BadIdx;
    VDJdata(UpdateIdx,:) = buildRefSeq(VDJdata(UpdateIdx,:),VDJheader,'germline','single');
    VDJdata(UpdateIdx,:) = updateVDJdata(VDJdata(UpdateIdx,:),VDJheader,Vmap,Dmap,Jmap);    
end

if nargout >=2
    varargout{1} = BadIdx;
end
