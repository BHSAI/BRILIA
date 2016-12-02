%trimGeneEdge will trim the gene ends according to several logic rules
%
% VDJdata = trimGeneEdge(VDJdata,NewHeader)
%
% 1) Identify "consensus" mismatch using majority rule, inclusive of 50%
% split.
%
% 2) Always trim at >= 3 consec mismatch and happens between 104C and 118W
% region.
%     EX: |||xxx||  --> oooooo|| if we are trimming left side, or
%                   --> |||ooooo if we are trimming right side.
% 
% 3) Always cut on concensus mismatches with the same NT for group >= 3
%     EX: 'ACGTGgGT'   will become   'ACGTG ggt'
%         'ACcTGgGT'   will become   'ACcTG ggt'
%         'AcGTGgGT'   will become   'AcGTG ggt'
%         'ACaTGgGT'   will become   'ACaTG ggt'
%         'ACGaGgGT'   will become   'ACGaG ggt'
%
% 4) Decide to trim for <= 2 consec mismatch based on how well it
% improves/worsen alignment score, for majority mismatch rule.
%     EX: |x||||||  --> oo|||||| 
%     EX: ||x|||||  --> ||x||||| no trim
%     EX: ||xx||||  --> oooo||||
%     EX: |x|x||||  --> oooo||||
%     EX: |||xx|||  --> |||xx||| no trim

function VDJdata = trimGeneEdge(VDJdata,NewHeader,varargin)
%Extract the VDJ database
if length(varargin) == 3
    Vmap = varargin{1};
    Dmap = varargin{2};
    Jmap = varargin{3};
else
    [Vmap, Dmap, Jmap] = getCurrentDatabase;
end
getHeaderVar;

%Look for better N
GrpNum = cell2mat(VDJdata(:,GrpNumLoc));
GrpNumUnq = unique(GrpNum);
IdxMap = 1:size(GrpNum,1);
for y = 1:length(GrpNumUnq)
    %Extract the group
    GrpLoc = (GrpNumUnq(y) == GrpNum);
    IdxLoc = IdxMap(GrpLoc);        
    Tdata = VDJdata(GrpLoc,:);
    
    %Makesure this if a functional sequence
    if strcmpi(Tdata{1,FunctLoc},'N')
        continue
    end
    
    %Find the location of consensus mismatch using the Classifier info
    Classifier = char(Tdata(:,FormClassLoc(2)));
    ConsMissCt = zeros(1,size(Classifier,2));
    for k = 1:size(Classifier,1)
        ConsMissCt = ConsMissCt + isstrprop(Classifier(k,:),'lower');
    end
    ConsMiss = ConsMissCt >= ceil(0.95*size(Classifier,1)); 
    
    %======================================================================
    %Extract the segment informations within the CDR3 regions
    Seq = Tdata{1,SeqLoc};
    VMDNJ = cell2mat(Tdata(1,LengthLoc));
    C104 = Tdata{1,CDR3Loc(3)};
    if C104 < 1
        C104 = 1;
    end
    W118 = Tdata{1,CDR3Loc(4)}; %This can be greater than SeqLength due to anchor point of J going over
    if W118 > length(Seq)
        W118 = length(Seq);
    end
    
    VdelCur = Tdata{1,DelLoc(1)};
    JdelCur = Tdata{1,DelLoc(end)};
    D5delCur = Tdata{1,DelLoc(2)};
    D3delCur = Tdata{1,DelLoc(3)};
    
    %Nucleotides at the VDJ edges
    Vnt = Seq(C104:VMDNJ(1));
    Mnt = Seq(VMDNJ(1)+1:sum(VMDNJ(1:2)));
    Dnt = Seq(sum(VMDNJ(1:2))+1:sum(VMDNJ(1:3)));
    MidLoc = floor(length(Dnt)/2);
    D5nt = Dnt(1:MidLoc);
    D3nt = Dnt(end-MidLoc+1:end);
    Nnt = Seq(sum(VMDNJ(1:3))+1:sum(VMDNJ(1:4)));
    Jnt = Seq(sum(VMDNJ(1:4))+1:W118);
    
    %Consens mismatched positions in the VDJ edges
    Vmiss = ConsMiss(C104:VMDNJ(1));
    Jmiss = ConsMiss(sum(VMDNJ(1:4))+1:W118);  
    Dmisses = ConsMiss(sum(VMDNJ(1:2))+1:sum(VMDNJ(1:3)));
    DmissL = Dmisses(1:MidLoc);
    DmissR = Dmisses(end-MidLoc+1:end);
    
    %For the Nvd regions, split based on p-nucleotide. EX: bmmmb = 12223.
    Mclass = Classifier(1,VMDNJ(1)+1:sum(VMDNJ(1:2)));
    Mloc = find(Mclass == 'm');
    if isempty(Mloc); Mloc = 0; end
    BMB = ones(1,length(Mclass))*2;
    BMB(1:Mloc(1)-1) = 1;
    BMB(Mloc(end)+1:end) = 3;
    
    %For the Ndj regions, split based on p-nucleotide. EX: pnnnp = 12223.
    Nclass = Classifier(1,sum(VMDNJ(1:3))+1:sum(VMDNJ(1:4)));
    Nloc = find(Nclass == 'n');
    if isempty(Nloc); Nloc = 0; end
    PNP = ones(1,length(Nclass))*2;
    PNP(1:Nloc(1)-1) = 1;
    PNP(Nloc(end)+1:end) = 3;
    
    %======================================================================
    %AUTOTRIM rule if you have consec mismatches within CDR3 region VDJs    
    MinConsec = 3; %Trim if you find this many consec mismatches    
    VtrimLen1 = findConsecLoc(Vmiss,MinConsec,'right');
    JtrimLen1 = findConsecLoc(Jmiss,MinConsec,'left');    
    D5trimLen1 = findConsecLoc(DmissL,MinConsec,'left');
    D3trimLen1 = findConsecLoc(DmissR,MinConsec,'right');
   
    %======================================================================
    %Trim edges if TDTscore is higher for trim current than its respective
    %Ref TDT.
    if sum(BMB==2) == 0
        Mscore = 0.5;
    else
        Mscore = calcTDTscore(Mnt(BMB==2));
    end
    
    if sum(PNP==2) == 0
        Nscore = 0.5;
    else
        Nscore = calcTDTscore(Nnt(PNP==2));
    end
    
    %Start with V edges
    Vloc = sort(find(Vmiss == 1),'descend');
    Vscores = zeros(length(Vloc),2);
    for v = 1:size(Vscores,1)
        Vscores(v,1) = Vloc(v);
        Vscores(v,2) = calcTDTscore([Vnt(Vloc(v):end) Mnt(BMB<=2)]);
    end
    ValidLoc = find(Vscores(:,2) > Mscore);
    if isempty(ValidLoc)
        VtrimLen2 = 0;
    else
        VtrimLen2 = length(Vnt) - Vscores(ValidLoc(end),1) + 1;
    end
    
    %Start with J edges
    Jloc = find(Jmiss == 1);
    Jscores = zeros(length(Jloc),2);
    for v = 1:size(Jscores,1)
        Jscores(v,1) = Jloc(v);
        Jscores(v,2) = calcTDTscore([Nnt(PNP>=2) Jnt(1:Jloc(v))]);
    end
    ValidLoc = find(Jscores(:,2) > Nscore);
    if isempty(ValidLoc)
        JtrimLen2 = 0;
    else
        JtrimLen2 = Jscores(ValidLoc(end),1);
    end
    
    %Start with D5 edges
    D5loc = find(DmissL == 1);
    D5scores = zeros(length(D5loc),2);
    for v = 1:size(D5scores,1)
        D5scores(v,1) = D5loc(v);
        D5scores(v,2) = calcTDTscore([Mnt(BMB>=2) D5nt(1:D5loc(v))]);
    end
    ValidLoc = find(D5scores(:,2) > Mscore);
    if isempty(ValidLoc)
        D5trimLen2 = 0;
    else
        D5trimLen2 = D5scores(ValidLoc(end),1);
    end
    
    %Start with D3 edges 
    D3loc = sort(find(DmissR == 1),'descend');
    D3scores = zeros(length(D3loc),2);
    for v = 1:size(D3scores,1)
        D3scores(v,1) = D3loc(v);
        D3scores(v,2) = calcTDTscore([D3nt(D3loc(v):end) Nnt(PNP<=2)]);
    end
    ValidLoc = find(D3scores(:,2) > Nscore);
    if isempty(ValidLoc)
        D3trimLen2 = 0;
    else
        D3trimLen2 = length(D3nt) - D3scores(ValidLoc(end),1) + 1;
    end
    
    %Cap the deletions by 5 nts...
    VtrimLens = [VtrimLen1 VtrimLen2];
    VtrimLens(VtrimLens > 5) = [];
    D5trimLens = [D5trimLen1 D5trimLen2];
    D5trimLens(D5trimLens > 5) = [];
    D3trimLens = [D3trimLen1 D3trimLen2];
    D3trimLens(D3trimLens > 5) = [];
    JtrimLens = [JtrimLen1 JtrimLen2];
    JtrimLens(JtrimLens > 5) = [];    
    
    %Determine the maximum cuts
    Vcut = max(VtrimLens); %max([VtrimLen1 VtrimLen2]);
    D5cut = max(D5trimLens); %max([D5trimLen1 D5trimLen2]);
    D3cut = max(D3trimLens); %max([D3trimLen1 D3trimLen2]);
    Jcut = max(JtrimLens); %max([JtrimLen1 JtrimLen2]);    
    
    %Correct only if there is a change
    if sum(Vcut + Jcut + D5cut + D3cut) > 0
        if VMDNJ(3)-D5cut-D3cut < 3; %Be cautious of the cut
            if D5trimLen1 == 0 %If there's no triple mismatch, ignore trim to preserve D
               D5cut = 0;
               disp('Prevented D5 cut');
            end
            
            if D3trimLen1 == 0 %If there's no triple mismatch, ignore trim to preserve D
                D3cut = 0;
                disp('Prevented D3 cut');
            end
            
            if VMDNJ(3)-D5cut-D3cut < 3 %If it's still too much, just ignore all trims
                disp('Prevented all D cut');
                D5cut = 0;
                D3cut = 0;
            end
        end

        VMDNJnew = [VMDNJ(1)-Vcut  VMDNJ(2)+Vcut+D5cut  VMDNJ(3)-D5cut-D3cut  VMDNJ(4)+Jcut+D3cut VMDNJ(5)-Jcut];
        NewDel = [VdelCur+Vcut  D5delCur+D5cut  D3delCur+D3cut  JdelCur+Jcut];

        if VMDNJnew(1)*VMDNJnew(3)*VMDNJnew(5) == 0; %Try not to lose the VDJ.
            continue; 
        end 
        if min(VMDNJnew) < 0; %Don't fix errors
            disp(['Weird D correction found in Group ' num2str(GrpNumUnq(y))]);
            continue; 
        end 
        
        Tdata(:,LengthLoc) = repmat(num2cell(VMDNJnew),size(Tdata,1),1);
        Tdata(:,DelLoc) = repmat(num2cell(NewDel),size(Tdata,1),1);    

        %Redo the alignments and classifier
        Tdata = buildVDJalignment(Tdata,NewHeader,Vmap,Dmap,Jmap);
        Tdata = buildRefSeq(Tdata,NewHeader,'same','germline','first');%Same length, germline substitution, on first sequence of each group
        Tdata = makeClassifier(Tdata,NewHeader);
        Tdata = appendMutCt(Tdata,NewHeader); %SHM infor on the VMDNJ segments
        
        VDJdata(IdxLoc,:) = Tdata;
    end
end    
