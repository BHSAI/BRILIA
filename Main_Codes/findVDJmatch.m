%findVDJmatch will look for the best V,D,J gene match for a sequence in the
%following order: find V preserving nts for DJ, find J perserving nts for
%D, find D. This is good for getting an initial guess VDJ annotation.
%
%  VDJdata = findVDJmatch(VDJdata,NewHeader)
%
%  VDJdata = findVDJmatch(VDJdata,NewHeader,Vmap,Dmap,Jmap)
%
%  VDJdata = findVDJmatch(VDJdata,NewHeader,Vmap,Dmap,Jmap,'update')
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
%    faster alignment
% 
%    After finding V, it will check to see if a 'TGG' of the 118W exists
%    after the 'TGT' of the conserved 104C, that is in frame. If it does
%    find one, it will force a J alignment such that the 118W are aligned.
%    This helps if the input seq is too short and ends exactly on the
%    conserved 118W. Otherwise, it either uses the seed from
%    seedCDR3positoin or a full J alignment to determine the J segment.
%
function [VDJdata,varargout] = findVDJmatch(VDJdata,NewHeader,varargin)
%Look if update was specified first
UpdateLoc = findCell(varargin,'update','MatchCase','any');
if UpdateLoc(1) ~= 0
    Update = 'update';
    varargin(UpdateLoc) = [];
else
    Update = '';
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
Dreserve = P.Results.Dreserve; %Preserve first 3 nts after V for D.

%Extract the VDJ database, if none is specified.
if isempty(Vmap) || isempty(Dmap) || isempty(Jmap)
    [Vmap, Dmap, Jmap] = getCurrentDatabase('change');
    [Vmap, Dmap, Jmap] = filterRefGene(Vmap,Dmap,Jmap);
end

%Bring some getHeaderVar variables out, since parfor can't handle it.
SeqLoc = findCell(NewHeader,{'nucleotide','Seq'});
LengthLoc = findCell(NewHeader,{'Length_V','Length_Nvd','Length_D','Length_Ndj','Length_J'});
DelLoc = findCell(NewHeader,{'V_Deletion3','D_Deletion5','D_Deletion3','J_Deletion5'});
FamNumLoc = findCell(NewHeader,{'V_MapNum','D_MapNum','J_MapNum'});
FamLoc = findCell(NewHeader,{'V_GeneName','D_GeneName','J_GeneName'});
MiscLoc = findCell(NewHeader,{'Misc'});
CDR3startLoc = findCell(NewHeader,{'CDR3_Start'});
CDR3endLoc = findCell(NewHeader,{'CDR3_End'});

%If given just a single sequence, then create a VDJmatrix.
if ~iscell(VDJdata)  %Process a single sequence, create a single-line VDJdata
    Seq = VDJdata;
    
    %Create the VDJdata default matrix
    [~, ~, HeaderData] = xlsread('Headers_BRILIA.xlsx'); %Obtain the VDJdata header info for output format
    NewHeader = HeaderData(2:end,1)';
    getHeaderVar;
    VDJdata = cell(size(Seq,1),length(NewHeader)); 
    VDJdata(:,TemplateLoc) = num2cell(ones(size(Seq,1),1)); %Always initialize TempCt column with 1.
    VDJdata(:,SeqNumLoc) = num2cell(1:size(Seq,1)); %Always assign a unique numbering order for VDJdata
    VDJdata{:,GrpNumLoc} = 1;
    VDJdata{:,SeqNameLoc} = '1';
    VDJdata{:,SeqLoc} = Seq;
end
    
BadIdx = zeros(size(VDJdata,1),1,'logical');
parfor j = 1:size(VDJdata,1)
    try
        Tdata = VDJdata(j,:);
        Seq = Tdata{1,SeqLoc};
        CDR3start = Tdata{1,CDR3startLoc};
        CDR3end = Tdata{1,CDR3endLoc};

        MissRate = 0.15; %Start with 15% miss rate for the V segment. Place inside parfor loop to reset.

        %Look for V gene for each seq
        Vnt = Seq(1:end-DJreserve); %Preserve last NTs to prevent matching V without D and J.
        AllowedMiss = ceil(MissRate * (length(Vnt)));
        Vmatch = findGeneMatch(Vnt,Vmap,'V',AllowedMiss,CDR3start);
        Vlen = sum(Vmatch{4}(1:2));
        MissRate = (Vlen - Vmatch{5}(1,1))/Vlen; %Re-establish MissRate for D and J
        
        %Look for J gene for each seq
        Jnt = Seq(Vlen+Dreserve+1:end); %Preserve D NTs for later matching, and pad seq for OverhangMatch success.
        AllowedMiss = ceil(MissRate * length(Jnt));  %Remember to shift CDR3end values
        CDR3endTemp = CDR3end - Vlen - Dreserve;
        CDR3endTemp(CDR3endTemp <= 0) = [];
        CDR3endPos = regexp(Jnt,'TGG');
        if ~isempty(CDR3endPos)
            CDR3endTemp = unique([CDR3endTemp(:); CDR3endPos(:)+2]);
        end
        Jmatch = findGeneMatch(Jnt,Jmap,'J',AllowedMiss,CDR3endTemp);
        Jlen = sum(Jmatch{4}(2:3));

        %Look for D gene for each seq
        Dnt = Seq(Vlen+1:end-Jlen);
        AllowedMiss = ceil(MissRate * length(Dnt));
        Dmatch = findGeneMatch(Dnt,Dmap,'D',AllowedMiss);

        %Extract VMDNJ lengths
        Vlmr = cell2mat(Vmatch(1,4));
        Jlmr = cell2mat(Jmatch(1,4));
        Dlmr = cell2mat(Dmatch(1,4));
        VMDNJ = [sum(Vlmr(1,1:2),2) Dlmr sum(Jlmr(1,2:3),2)];
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

        %Perform quality control check
        
        VDJdata(j,:) = Tdata;
    catch
        BadIdx(j) = 1;
    end
end

%If there are errors, show them now.
BadLoc = find(BadIdx == 1);
for b = 1:length(BadLoc)
    ErrorMsg = sprintf('Errored at %s, sequence # %d',mfilename,BadLoc(b));
    disp(ErrorMsg);
    VDJdata{BadLoc(b),MiscLoc} = ErrorMsg;
end

%Update VDJdata 
if strcmpi(Update,'update')
    UpdateIdx = ~BadIdx;
    VDJdata(UpdateIdx,:) = buildRefSeq(VDJdata(UpdateIdx,:),NewHeader,'germline','single');
    VDJdata(UpdateIdx,:) = updateVDJdata(VDJdata(UpdateIdx,:),NewHeader,Vmap,Dmap,Jmap);    
end

if nargout >=2
    varargout{1} = BadIdx;
end