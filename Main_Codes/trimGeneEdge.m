%trimGeneEdge will trim the V-D-J gene ends according to several logic
%rules, which can be added or modified over time depending on new
%knowledge. This is the function to add exception rules to gene edge
%finding.
%
%  VDJdata = trimGeneEdge(VDJdata,NewHeader)
%
%  NOTES ON LOGIC RULES
%    1) Identify "consensus" mismatch using 50% majority rule, inclusive of
%       50% split.
%
%    2) Trim at >= 3 consec mismatches that happens between 104C and 118W.
%       |||xxx||  ->  oooooo|| if we are trimming left side.
%                 ->  |||ooooo if we are trimming right side.
% 
%    3) Trim at consensus mismatch if it increases TDTscore. For lack of initial N region, the TDTscore starts at 0.5.
%       'ACGTGgGT       'ACGTG ggt'
%       'ACcTGgGT'  ->  'ACcTG ggt'
%       'AcGTGgGT'      'AcGTG ggt'
%
%  See also BRILIA, findVDJmatch

function VDJdata = trimGeneEdge(VDJdata,NewHeader,varargin)
%Extract the VDJ database
if length(varargin) >= 3
    [Vmap,Dmap,Jmap] = deal(varargin{1:3});
else
    [Vmap,Dmap,Jmap] = getCurrentDatabase;
end
getHeaderVar;

%Look for better N
GrpNum = cell2mat(VDJdata(:,GrpNumLoc));
UnqGrpNum = unique(GrpNum);
UpdateIdx = zeros(size(VDJdata,1),1,'logical');
for y = 1:length(UnqGrpNum)
    try
        %Extract the basic info
        GrpIdx = find(UnqGrpNum(y) == GrpNum);
        Tdata = VDJdata(GrpIdx,:); 
        RefSeq = Tdata{1,RefSeqLoc};
        VMDNJ = cell2mat(Tdata(1,LengthLoc));
        [VdelCur,D5delCur,D3delCur,JdelCur] = deal(Tdata{1,DelLoc});
        [CDR3start,CDR3end] = deal(Tdata{1,CDR3Loc(3:4)});
        if CDR3start < 1
            CDR3start = 1;
        end
        if CDR3end > length(RefSeq)
            CDR3end = length(RefSeq);
        end
        
        %Extract necessary V informations
        VmapNum = Tdata{1,FamNumLoc(1)};
        VallowedDel = Vmap{VmapNum(1),end} - 3; %Correct -3 as deletion length is AFTER C codon.
        if VallowedDel < 0; VallowedDel = 25; end        
 
        %Find the mismatched nts with respect to 1st seq of cluster only.
        XlocRef = RefSeq == 'X';
        ConsMissCt = zeros(size(RefSeq));
        for k = 1:size(Tdata,1)
            Seq = Tdata{k,SeqLoc};
            XlocSeq = Seq == 'X';
            MissLoc = ~(RefSeq == Seq | XlocSeq | XlocRef);
            ConsMissCt = ConsMissCt + MissLoc;
        end    
        
        %Find consensus misses that are >= MaxMiss as seen in the V
        %segment, which is the most accurate segment
        VconsSeg = VMDNJ(1)+VdelCur-VallowedDel;
        MaxMiss = max(ConsMissCt(1:VconsSeg));
        if size(Tdata,1) == 1 || isempty(MaxMiss) || MaxMiss == 0
            ConsMiss = ConsMissCt >= 1; %Must be atleast 1
        else
            ConsMiss = ConsMissCt >= MaxMiss;
        end
               
        %Nucleotides at the VDJ edges
        Vnt = Seq(CDR3start:VMDNJ(1));
        Mnt = Seq(VMDNJ(1)+1:sum(VMDNJ(1:2)));
        Dnt = Seq(sum(VMDNJ(1:2))+1:sum(VMDNJ(1:3)));
        Nnt = Seq(sum(VMDNJ(1:3))+1:sum(VMDNJ(1:4)));
        Jnt = Seq(sum(VMDNJ(1:4))+1:CDR3end);

        %Consens mismatched positions in the VDJ edges
        Vmiss = ConsMiss(CDR3start:VMDNJ(1));
        Dmiss = ConsMiss(sum(VMDNJ(1:2))+1:sum(VMDNJ(1:3)));
        Jmiss = ConsMiss(sum(VMDNJ(1:4))+1:CDR3end);  

        %Split the Dnts and Dmiss into 5' and 3' side
        MidLoc = floor(length(Dnt)/2);
        Dnt5 = Dnt(1:MidLoc);
        Dnt3 = Dnt(end-MidLoc+1:end);
        Dmiss5 = Dmiss(1:MidLoc);
        Dmiss3 = Dmiss(end-MidLoc+1:end);
        
        %Determine P nt locations. Used to calculate correct Nscore.
        VsideP = findPnts(RefSeq,[1,VMDNJ(1)],'right',0,VdelCur);
        DsideP = findPnts(RefSeq,[sum(VMDNJ(1:2))+1,sum(VMDNJ(1:3))],'both',D5delCur,D3delCur);
        JsideP = findPnts(RefSeq,[sum(VMDNJ(1:4))+1,sum(VMDNJ)],'left',0,JdelCur);
        PntLoc = VsideP & DsideP & JsideP;
        
        %Track p and n nucleotides  using number scheme [1 1 2 2 2 3 3 3],
        %where 1 is left side p-nt, 2 is n-nt, and 3 is right side p-nt.
        PntMat = double(~PntLoc)*2;
        
        BMB = PntMat(VMDNJ(1)+1:sum(VMDNJ(1:2)));
        Bidx = find(BMB == 2);
        if ~isempty(Bidx)
            BMB(1:Bidx(1)-1) = 1;
            BMB(Bidx(end)+1:end) = 3;        
        end
        
        PNP = PntMat(sum(VMDNJ(1:3))+1:sum(VMDNJ(1:4)));
        Pidx = find(PNP == 2);
        if ~isempty(Pidx)
            PNP(1:Pidx(1)-1) = 1;
            PNP(Pidx(end)+1:end) = 3;
        end
        

        %======================================================================
        %AUTOTRIM rule if you have consec mismatches within CDR3 region VDJs    
        MinConsec = 3; %Trim if you find this many consec mismatches    
        VtrimLen1 = findConsecLoc(Vmiss,MinConsec,'right');
        JtrimLen1 = findConsecLoc(Jmiss,MinConsec,'left');    
        D5trimLen1 = findConsecLoc(Dmiss5,MinConsec,'left');
        D3trimLen1 = findConsecLoc(Dmiss3,MinConsec,'right');

        %======================================================================
        %Only trim edges that increase the TDT score
        
        %Calculate the baseline TDTscores
        if isempty(BMB)
            Mscore = 0.5;
        elseif sum(BMB==2) == 0
            Mscore = 0.5;
        else
            Mscore = calcTDTscore(Mnt(BMB==2));
        end

        if isempty(PNP)
            Nscore = 0.5;
        elseif max(PNP==2) == 0
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
        D5loc = find(Dmiss5 == 1);
        D5scores = zeros(length(D5loc),2);
        for v = 1:size(D5scores,1)
            D5scores(v,1) = D5loc(v);
            D5scores(v,2) = calcTDTscore([Mnt(BMB>=2) Dnt5(1:D5loc(v))]);
        end
        ValidLoc = find(D5scores(:,2) > Mscore);
        if isempty(ValidLoc)
            D5trimLen2 = 0;
        else
            D5trimLen2 = D5scores(ValidLoc(end),1);
        end

        %Start with D3 edges 
        D3loc = sort(find(Dmiss3 == 1),'descend');
        D3scores = zeros(length(D3loc),2);
        for v = 1:size(D3scores,1)
            D3scores(v,1) = D3loc(v);
            D3scores(v,2) = calcTDTscore([Dnt3(D3loc(v):end) Nnt(PNP<=2)]);
        end
        ValidLoc = find(D3scores(:,2) > Nscore);
        if isempty(ValidLoc)
            D3trimLen2 = 0;
        else
            D3trimLen2 = length(Dnt3) - D3scores(ValidLoc(end),1) + 1;
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
        Vcut = max(VtrimLens); 
        D5cut = max(D5trimLens);
        D3cut = max(D3trimLens);
        Jcut = max(JtrimLens); 

        %Correct only if there is a change
        if sum(Vcut + Jcut + D5cut + D3cut) > 0
            %Be cautious of trimming D's too much
            if VMDNJ(3)-D5cut-D3cut < 3;
                %If there's no triple mismatch, ignore trim to preserve D
                if D5trimLen1 == 0 
                    D5cut = 0;
                    %disp('Prevented D5 cut');
                end
                
                %If there's no triple mismatch, ignore trim to preserve D
                if D3trimLen1 == 0
                    D3cut = 0;
                    %disp('Prevented D3 cut');
                end
                
                %If it's still too much, just ignore all trims
                if VMDNJ(3)-D5cut-D3cut < 3 
                    %disp('Prevented all D cut');
                    D5cut = 0;
                    D3cut = 0;
                end
            end

            VMDNJnew = [VMDNJ(1)-Vcut  VMDNJ(2)+Vcut+D5cut  VMDNJ(3)-D5cut-D3cut  VMDNJ(4)+Jcut+D3cut VMDNJ(5)-Jcut];
            NewDel = [VdelCur+Vcut  D5delCur+D5cut  D3delCur+D3cut  JdelCur+Jcut];

            if VMDNJnew(1)*VMDNJnew(3)*VMDNJnew(5) == 0; %Try not to lose the VDJ.
               % continue; 
            end 
            if min(VMDNJnew) < 0; %Don't fix errors
                warning('Warning at %s, sequence group # %d, negative VMDNJ',mfilename,UnqGrpNum(y));
               % continue; 
            end

            %Update the necessary fields 
            Tdata(:,LengthLoc) = repmat(num2cell(VMDNJnew),size(Tdata,1),1);
            Tdata(:,DelLoc) = repmat(num2cell(NewDel),size(Tdata,1),1);    
            VDJdata(GrpIdx,:) = Tdata;
            UpdateIdx(GrpIdx) = 1;
        end
    catch
        warning('Warning at %s, sequence group # %d',mfilename,UnqGrpNum(y));
        pause        
    end
end

%Update those that have changed
VDJdata(UpdateIdx,:) = buildRefSeq(VDJdata(UpdateIdx,:),NewHeader,'germline','first'); %must do first seq of all cluster
VDJdata(UpdateIdx,:) = updateVDJdata(VDJdata(UpdateIdx,:),NewHeader,varargin);