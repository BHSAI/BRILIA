%trimGeneEdge will trim the V-(D)-J gene ends according to several logic
%rules, which can be added or modified over time depending on newly gained
%knowledge. This is the function to add exception rules to gene edge
%finding.
%
%  VDJdata = trimGeneEdge(VDJdata, VDJheader)
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
%  INPUT
%    VDJdata: main BRILIA data cell
%    VDJheader: main BRILIA header cell
%    DB: Gene database structure (getCurrentDatabase.m)
%
%  OUTPUT
%    VDJdata: modified VDJdata with more precise Nvd, Ndj, or Nvj.
%
%  See also BRILIA, findVDJmatch

function VDJdata = trimGeneEdge(VDJdata, Map, DB)
if isempty(VDJdata) ; return; end

%Look for better N
if ~iscell(VDJdata{1}) %VDJdata is not spliced
    GrpNum = cell2mat(VDJdata(:, Map.GrpNum));
    UnqGrpNum = unique(GrpNum);
    for y = 1:length(UnqGrpNum)
        if ~mod(y, 1000)
            showStatus(sprintf('  Trimming genes  %d / %d.', y, length(UnqGrpNum)));
        end
        IdxLoc = UnqGrpNum(y) == GrpNum;
        VDJdata(IdxLoc, :) = trimGeneEdgePerGroup(VDJdata(IdxLoc, :), Map, DB);
    end
else %VDJdata is spliced for parfor
    parfor y = 1:length(VDJdata)
        VDJdata{y} = trimGeneEdgePerGroup(VDJdata{y}, Map, DB);
    end
end

function Tdata = trimGeneEdgePerGroup(Tdata, Map, DB)
UpdateThis = 0;
if contains(Map.Chain, 'H', 'ignorecase', true)
    %Extract the basic info
    RefSeq = Tdata{1, Map.hRefSeq};
    VMDNJ = cell2mat(Tdata(1, Map.hLength));
    [VdelCur, D5delCur, D3delCur, JdelCur] = deal(Tdata{1, Map.hDel});
    [CDR3s, CDR3e] = deal(Tdata{1, Map.hCDR3(3:4)});

    %Ensure all info is available before trimming
    if ~(isempty(RefSeq) || isempty(VMDNJ) || isempty(VdelCur) || isempty(D5delCur) || isempty(D3delCur) || isempty(JdelCur) || isempty(CDR3s) || isempty(CDR3e))
        %Ensure CDR3s and e span the range of RefSeq
        if CDR3s < 1
            CDR3s = 1;
        elseif CDR3s+2 > length(RefSeq)
            warning('%s: start of heavy cdr3 is beyond the sequence\n  Seq = %s\n  RefSeq = %s\n  CDR3s = %d\n  CDR3e = %d\n', mfilename, Tdata{1, Map.hSeq}, RefSeq, CDR3s, CDR3e);
            return
        end
        if CDR3e > length(RefSeq)
            CDR3e = length(RefSeq);
        end

        %Find the mismatched nts with respect to 1st seq of cluster only.
        XlocRef = RefSeq == 'N';
        ConsMissCt = zeros(size(RefSeq));
        for k = 1:size(Tdata, 1)
            Seq = Tdata{k, Map.hSeq};
            if length(Seq) ~= length(RefSeq) 
                return
            end
            if isempty(Seq); return; end
            XlocSeq = Seq == 'N';
            MissLoc = ~(RefSeq == Seq | XlocSeq | XlocRef);
            ConsMissCt = ConsMissCt + MissLoc;
        end

        %Find consensus misses that are >= MaxMiss as seen in the V
        %segment, which is the most accurate segment
        MaxMiss = max(ConsMissCt(1:CDR3s+2));
        if size(Tdata, 1) == 1 || isempty(MaxMiss) || MaxMiss == 0
            ConsMiss = ConsMissCt >= 1; %Must be atleast 1
        else
            ConsMiss = ConsMissCt >= MaxMiss;
        end

        %Nucleotides at the VDJ edges
        Vnt = Seq(CDR3s:VMDNJ(1));
        Mnt = Seq(VMDNJ(1)+1:sum(VMDNJ(1:2)));
        Dnt = Seq(sum(VMDNJ(1:2))+1:sum(VMDNJ(1:3)));
        Nnt = Seq(sum(VMDNJ(1:3))+1:sum(VMDNJ(1:4)));
        Jnt = Seq(sum(VMDNJ(1:4))+1:CDR3e);

        %Consens mismatched positions in the VDJ edges
        Vmiss = ConsMiss(CDR3s:VMDNJ(1));
        Dmiss = ConsMiss(sum(VMDNJ(1:2))+1:sum(VMDNJ(1:3)));
        Jmiss = ConsMiss(sum(VMDNJ(1:4))+1:CDR3e);  

        %Split the Dnts and Dmiss into 5' and 3' side
        MidLoc = floor(length(Dnt)/2);
        Dnt5 = Dnt(1:MidLoc);
        Dnt3 = Dnt(end-MidLoc+1:end);
        Dmiss5 = Dmiss(1:MidLoc);
        Dmiss3 = Dmiss(end-MidLoc+1:end);

        %Determine P nt locations. Used to calculate correct Nscore.
        VsideP = findPnts(RefSeq, [1, VMDNJ(1)], 'right', 0, VdelCur);
        DsideP = findPnts(RefSeq, [sum(VMDNJ(1:2))+1, sum(VMDNJ(1:3))], 'both', D5delCur, D3delCur);
        JsideP = findPnts(RefSeq, [sum(VMDNJ(1:4))+1, sum(VMDNJ)], 'left', 0, JdelCur);
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
        VtrimLen1 = findConsecLoc(Vmiss, MinConsec, 'right');
        JtrimLen1 = findConsecLoc(Jmiss, MinConsec, 'left');    
        D5trimLen1 = findConsecLoc(Dmiss5, MinConsec, 'left');
        D3trimLen1 = findConsecLoc(Dmiss3, MinConsec, 'right');

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
        Vloc = sort(find(Vmiss == 1), 'descend');
        Vscores = zeros(length(Vloc), 2);
        for v = 1:size(Vscores, 1)
            Vscores(v, 1) = Vloc(v);
            Vscores(v, 2) = calcTDTscore([Vnt(Vloc(v):end) Mnt(BMB<=2)]);
        end
        ValidLoc = find(Vscores(:, 2) > Mscore);
        if isempty(ValidLoc)
            VtrimLen2 = 0;
        else
            VtrimLen2 = length(Vnt) - Vscores(ValidLoc(end), 1) + 1;
        end

        %Start with J edges
        Jloc = find(Jmiss == 1);
        Jscores = zeros(length(Jloc), 2);
        for v = 1:size(Jscores, 1)
            Jscores(v, 1) = Jloc(v);
            Jscores(v, 2) = calcTDTscore([Nnt(PNP>=2) Jnt(1:Jloc(v))]);
        end
        ValidLoc = find(Jscores(:, 2) > Nscore);
        if isempty(ValidLoc)
            JtrimLen2 = 0;
        else
            JtrimLen2 = Jscores(ValidLoc(end), 1);
        end

        %Start with D5 edges
        D5loc = find(Dmiss5 == 1);
        D5scores = zeros(length(D5loc), 2);
        for v = 1:size(D5scores, 1)
            D5scores(v, 1) = D5loc(v);
            D5scores(v, 2) = calcTDTscore([Mnt(BMB>=2) Dnt5(1:D5loc(v))]);
        end
        ValidLoc = find(D5scores(:, 2) > Mscore);
        if isempty(ValidLoc)
            D5trimLen2 = 0;
        else
            D5trimLen2 = D5scores(ValidLoc(end), 1);
        end

        %Start with D3 edges 
        D3loc = sort(find(Dmiss3 == 1), 'descend');
        D3scores = zeros(length(D3loc), 2);
        for v = 1:size(D3scores, 1)
            D3scores(v, 1) = D3loc(v);
            D3scores(v, 2) = calcTDTscore([Dnt3(D3loc(v):end) Nnt(PNP<=2)]);
        end
        ValidLoc = find(D3scores(:, 2) > Nscore);
        if isempty(ValidLoc)
            D3trimLen2 = 0;
        else
            D3trimLen2 = length(Dnt3) - D3scores(ValidLoc(end), 1) + 1;
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
            if VMDNJ(3)-D5cut-D3cut < 3
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

            NewVMDNJ = [VMDNJ(1)-Vcut  VMDNJ(2)+Vcut+D5cut  VMDNJ(3)-D5cut-D3cut  VMDNJ(4)+Jcut+D3cut VMDNJ(5)-Jcut];
            NewDel = [VdelCur+Vcut  D5delCur+D5cut  D3delCur+D3cut  JdelCur+Jcut];

            %Check to ensure the VMDNJ and del makes sense
            if NewVMDNJ(1)*NewVMDNJ(3)*NewVMDNJ(5) == 0 %Try not to lose the VDJ.
               return 
            end 
            if min(NewVMDNJ) < 0 %Don't fix errors
               warning('%s: Sequence group # %d has negative VMDNJ', mfilename, UnqGrpNum(y));
               return 
            end

            %Update the necessary fields 
            Tdata(:, Map.hLength) = repmat(num2cell(NewVMDNJ), size(Tdata, 1), 1);
            Tdata(:, Map.hDel) = repmat(num2cell(NewDel), size(Tdata, 1), 1);
            UpdateThis = 1;
        end
    end %End of valid info 
end %End of heavy chain fix

if contains(Map.Chain, 'L', 'ignorecase', true)     
    %Extract the basic info
    RefSeq = Tdata{1, Map.lRefSeq};
    VNJ = cell2mat(Tdata(1, Map.lLength));
    [VdelCur, JdelCur] = deal(Tdata{1, Map.lDel});
    [CDR3s, CDR3e] = deal(Tdata{1, Map.lCDR3(3:4)});

    %Ensure all info is available before trimming
    if ~(isempty(RefSeq) || isempty(VNJ) || isempty(VdelCur) || isempty(JdelCur) || isempty(CDR3s) || isempty(CDR3e))
    
        %Ensure CDR3s and e span the range of RefSeq
        if CDR3s < 1
            CDR3s = 1;
        elseif CDR3s+2 > length(RefSeq)
            fprintf('WARNING:%s: start of light cdr3 is beyond the sequence\n  Seq = %s\n  RefSeq = %s\n  CDR3s = %d\n  CDR3e = %d\n', mfilename, Tdata{1, Map.lSeq}, RefSeq, CDR3s, CDR3e);
            return
        end
        if CDR3e > length(RefSeq)
            CDR3e = length(RefSeq);
        end

        %Find the mismatched nts with respect to 1st seq of cluster only.
        XlocRef = RefSeq == 'N';
        ConsMissCt = zeros(size(RefSeq));
        for k = 1:size(Tdata, 1)
            Seq = Tdata{k, Map.lSeq};
            if isempty(Seq); return; end
            XlocSeq = Seq == 'N';
            MissLoc = ~(RefSeq == Seq | XlocSeq | XlocRef);
            ConsMissCt = ConsMissCt + MissLoc;
        end

        %Find consensus misses that are >= MaxMiss as seen in the V
        %segment, which is the most accurate segment
        MaxMiss = max(ConsMissCt(1:CDR3s+2));
        if size(Tdata, 1) == 1 || isempty(MaxMiss) || MaxMiss == 0
            ConsMiss = ConsMissCt >= 1; %Must be atleast 1
        else
            ConsMiss = ConsMissCt >= MaxMiss;
        end

        %Nucleotides at the VDJ edges
        Vnt = Seq(CDR3s:VNJ(1));
        Nnt = Seq(VNJ(1)+1:sum(VNJ(1:2)));
        Jnt = Seq(sum(VNJ(1:2))+1:CDR3e);

        %Consens mismatched positions in the VDJ edges
        Vmiss = ConsMiss(CDR3s:VNJ(1));
        Jmiss = ConsMiss(sum(VNJ(1:2))+1:CDR3e);  

        %Determine P nt locations. Used to calculate correct Nscore.
        VsideP = findPnts(RefSeq, [1, VNJ(1)], 'right', 0, VdelCur);
        JsideP = findPnts(RefSeq, [sum(VNJ(1:2))+1, sum(VNJ)], 'left', 0, JdelCur);
        PntLoc = VsideP & JsideP;

        %Track p and n nucleotides  using number scheme [1 1 2 2 2 3 3 3], 
        %where 1 is left side p-nt, 2 is n-nt, and 3 is right side p-nt.
        PntMat = double(~PntLoc)*2;
        PNP = PntMat(VNJ(1)+1:sum(VNJ(1:2)));
        Pidx = find(PNP == 2);
        if ~isempty(Pidx)
            PNP(1:Pidx(1)-1) = 1;
            PNP(Pidx(end)+1:end) = 3;
        end

        %======================================================================
        %AUTOTRIM rule if you have consec mismatches within CDR3 region VDJs    
        MinConsec = 3; %Trim if you find this many consec mismatches    
        VtrimLen1 = findConsecLoc(Vmiss, MinConsec, 'right');
        JtrimLen1 = findConsecLoc(Jmiss, MinConsec, 'left');

        %======================================================================
        %Only trim edges that increase the TDT score

        %Calculate the baseline TDTscores
        if isempty(PNP)
            Nscore = 0.5;
        elseif max(PNP==2) == 0
            Nscore = 0.5;
        else
            Nscore = calcTDTscore(Nnt(PNP==2));
        end

        %Start with V edges
        Vloc = sort(find(Vmiss == 1), 'descend');
        Vscores = zeros(length(Vloc), 2);
        for v = 1:size(Vscores, 1)
            Vscores(v, 1) = Vloc(v);
            Vscores(v, 2) = calcTDTscore([Vnt(Vloc(v):end) Nnt(PNP<=2)]);
        end
        ValidLoc = find(Vscores(:, 2) > Nscore);
        if isempty(ValidLoc)
            VtrimLen2 = 0;
        else
            VtrimLen2 = length(Vnt) - Vscores(ValidLoc(end), 1) + 1;
        end

        %Start with J edges
        Jloc = find(Jmiss == 1);
        Jscores = zeros(length(Jloc), 2);
        for v = 1:size(Jscores, 1)
            Jscores(v, 1) = Jloc(v);
            Jscores(v, 2) = calcTDTscore([Nnt(PNP>=2) Jnt(1:Jloc(v))]);
        end
        ValidLoc = find(Jscores(:, 2) > Nscore);
        if isempty(ValidLoc)
            JtrimLen2 = 0;
        else
            JtrimLen2 = Jscores(ValidLoc(end), 1);
        end

        %Cap the deletions by 5 nts...
        VtrimLens = [VtrimLen1 VtrimLen2];
        VtrimLens(VtrimLens > 5) = [];
        JtrimLens = [JtrimLen1 JtrimLen2];
        JtrimLens(JtrimLens > 5) = [];    

        %Determine the maximum cuts
        Vcut = max(VtrimLens); 
        Jcut = max(JtrimLens); 

        %Correct only if there is a change
        if sum(Vcut + Jcut) > 0
            NewVNJ = [VNJ(1)-Vcut  VNJ(2)+Vcut+Jcut VNJ(3)-Jcut];
            NewDel = [VdelCur+Vcut  JdelCur+Jcut];

            %Check to ensure the VMDNJ and del makes sense
            if NewVNJ(1)*NewVNJ(3) == 0 %Try not to lose the VDJ.
               return 
            end 
            if min(NewVNJ) < 0 %Don't fix errors
               warning('%s: Sequence group # %d has negative VNJ', mfilename, UnqGrpNum(y));
               return 
            end

            %Update the necessary fields 
            Tdata(:, Map.lLength) = repmat(num2cell(NewVNJ), size(Tdata, 1), 1);
            Tdata(:, Map.lDel) = repmat(num2cell(NewDel), size(Tdata, 1), 1);    
            UpdateThis = 1;
        end
    end %End of valid info 
end %End of Light chain fix    

if UpdateThis
    Tdata = buildRefSeq(Tdata, Map, DB, Map.Chain, 'germline', 'first'); %must do first seq of all cluster
    Tdata = updateVDJdata(Tdata, Map, DB);
end