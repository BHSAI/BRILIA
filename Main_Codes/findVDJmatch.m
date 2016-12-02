%findVDJmatch will look for the best V,D,J gene match in a sequence, in the
%following order: find V preserving nts for DJ, find J perserving nts for
%D, find D. This is good for getting an initial guess, or for redoing a
%match after indel correction.

function VDJdata = findVDJmatch(VDJdata,NewHeader,varargin)
PoolName = gcp('nocreate');
if size(VDJdata,1) > 200 && isempty(PoolName) == 1
    CoreNum = feature('numCores');
    PoolName = parpool(CoreNum,'AttachedFiles',{'getHeaderVar.m'});
end

%Extract the VDJ database
if length(varargin) == 3
    Vmap = varargin{1};
    Dmap = varargin{2};
    Jmap = varargin{3};
else
    [Vmap, Dmap, Jmap] = getCurrentDatabase;
end
getHeaderVar;

if ~iscell(VDJdata)  %Process a single sequence, create a single-line VDJdata
    Seq = VDJdata;
    %Extract relevant column locations
    VDJdata = cell(1,length(NewHeader));
    VDJdata{1,SeqNumLoc} = 1;
    VDJdata{1,GrpNumLoc} = 1;
    VDJdata{1,SeqLoc} = Seq;
end

%Set number of nts to preserve when match V then J the D.
DJreserve = 9; %Preserve last 9 for J.
Dreserve = 3;  %Preserve first 3 nts after V for D.

%Determine V and J pad lengths, based on allowed deletion counts. Helps
%with overhang matches by convolveSeq.
VpadCt = max(cell2mat(Vmap(:,end))) - DJreserve;
JpadCt = max(cell2mat(Jmap(:,end))) - Dreserve;
Vpad = repmat('*',1,VpadCt);
Jpad = repmat('*',1,JpadCt);

%Bring some getHeaderVar variables out, since parfor can't handle it.
SeqLoc = findHeader(NewHeader,{'nucleotide','Seq'});
    SeqLoc(SeqLoc==0) = []; SeqLoc = SeqLoc(1); 
LengthLoc = findHeader(NewHeader,{'vAlignLength' 'n2AlignLength' 'dAlignLength' 'n1AlignLength' 'jAlignLength','Length_V','Length_Nvd','Length_D','Length_Ndj','Length_J'});
    LengthLoc(LengthLoc==0) = [];
DelLoc = findHeader(NewHeader,{'vDeletion','d5Deletion','d3Deletion','jDeletion','V_Deletion3','D_Deletion5','D_Deletion3','J_Deletion5'});
    DelLoc(DelLoc==0) = [];
FamNumLoc = findHeader(NewHeader,{'vMapNum','dMapNum','jMapNum','V_MapNum','D_MapNum','J_MapNum'});
    FamNumLoc(FamNumLoc==0) = [];
FamLoc = findHeader(NewHeader,{'vMaxResolved','dMaxResolved','jMaxResolved','V_GeneName','D_GeneName','J_GeneName'});
    FamLoc(FamLoc==0) = []; 

parfor j = 1:size(VDJdata,1)
    Tdata = VDJdata(j,:);
    
    Seq = Tdata{1,SeqLoc};
    MissRate = 0.15; %Start with 15% miss rate for the V segment
    
    %Look for V gene for each seq
    Vnt = [Seq(1:end-DJreserve) Vpad]; %Preserve last NTs to prevent matching V without D and J.
    AllowedMiss = ceil(MissRate * (length(Vnt) - VpadCt));   
    Vmatch = findGeneMatch(Vnt,Vmap,'V',AllowedMiss); %Allow full leniency to work.
    Vmatch{4}(1,3) = Vmatch{4}(1,3) - length(Vpad); %Correct sample LMR for the padded NT match        
    Vlen = sum(Vmatch{4}(1:2));
    
    MissRate = (Vlen - Vmatch{5}(1,1))/Vlen; %Re-establish MissRate for D and J
    
    %Look for J gene for each seq
    Jnt = [Jpad Seq(Vlen+Dreserve+1:end)]; %Preserve D NTs for later matching, and pad seq for OverhangMatch success.
    AllowedMiss = ceil(MissRate * (length(Jnt) - JpadCt));
    Jmatch = findGeneMatch(Jnt,Jmap,'J',AllowedMiss);
    Jmatch{4}(1,1) = Jmatch{4}(1,1) - length(Jpad); %Correct sample LMR for the padded NT match
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
    
    VDJdata(j,:) = Tdata;
end

%Fill in the details now
VDJdata = buildRefSeq(VDJdata,NewHeader,'germline','single'); %must do singles, since group therapy not done.
VDJdata = buildVDJalignment(VDJdata,NewHeader,Vmap,Dmap,Jmap); %Alignment Info
VDJdata = makeClassifier(VDJdata,NewHeader); %Classifier + FormattedSeq
VDJdata = appendMutCt(VDJdata,NewHeader); %SHM infor on the VMDNJ segments
VDJdata = findCDR3(VDJdata,NewHeader); %Get the CDR3 seq and info 