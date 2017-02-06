%alignVbyC will align all V genes by the conserved C codon. Pseudogenes are
%considered, no C's are not . Returns all V, Vnumber on Vmap, and Vname, in
%aligned formed. Good for using this to reduce V search speed by looking
%for the V gene up to C, and then of those that matched, find the exact V
%gene match.

function Vmap = findVgeneC(Vmap)
%Get the functional genes
ValidV = zeros(size(Vmap,1),1) > 1;
for j = 1:size(ValidV,1);
    if strcmpi(Vmap{j,7},'f')
        ValidV(j) = 1;
    end
end
VmapT = Vmap(ValidV,:);

%Segregate by gene family numbers
[UnqFamNum,~,UnqIdx] = unique(VmapT(:,4));
ConservedSeq = cell(size(UnqFamNum,1),1);
for y = 1:length(UnqFamNum)
    UnqIdxLoc = find(UnqIdx == y);   
    RefSeq = VmapT{UnqIdxLoc(1),1};
    AlignMap = zeros(length(UnqIdxLoc),2);
    MaxLen = 0;
    for j = 1:length(UnqIdxLoc)
        CurSeq = VmapT{UnqIdxLoc(j),1};
        [~, SDF, StartAt, ~] = convolveSeq(CurSeq,RefSeq,0,0);
        AlignMap(j,:) = StartAt';
        if length(VmapT{UnqIdxLoc(j),1}) > MaxLen
            MaxLen = length(VmapT{UnqIdxLoc(j),1});
        end
    end
    
    %Creating alignment text
    MaxPad = abs(min(AlignMap(AlignMap(:,2)<0,2)));
    if isempty(MaxPad)
        MaxPad = 0;
    end
    AlignTxt = repmat(repmat('X',1,MaxLen+MaxPad),length(UnqIdxLoc),1);
    
    for j = 1:length(UnqIdxLoc)
        CurSeq = VmapT{UnqIdxLoc(j),1};
        
        StartAt = AlignMap(j,2) + MaxPad;
        if StartAt <= 0;
            StartAt = abs(StartAt) + 1;
        end
        EndAt = StartAt + length(CurSeq) - 1;

        AlignTxt(j,StartAt:EndAt) = CurSeq;
    end
    
    %Creating the alignment consesnsus count
    BctNorm = seqprofile(AlignTxt,'ALPHABET','NT','AMBIGUOUS','IGNORE','COUNT','FALSE');
    BctCons95 = max(BctNorm,[],1)>=0.50;
    MaxNT = zeros(1,size(BctNorm,2));
    for j = 1:length(MaxNT)
        CurMax = find(BctNorm(:,j) == max(BctNorm(:,j)));
        MaxNT(j) = CurMax(1);
    end
    ConsSeq = int2nt(MaxNT);
    
    %Convert to AA, find best frame
    ConsAA = nt2aa(ConsSeq,'frame','all','acgtonly','false');
    BestRF = zeros(3,2);
    for f = 1:3
        BestRF(f,2) = sum(ConsAA{f} == '*');
        BestRF(f,1) = f;
    end
    BestRF = sortrows(BestRF,2);
    RF = BestRF(1,1);
           
    %Find the conserved C "TGT" coding
    Clocs = find(ConsAA{RF} == 'C');
    Clocs = Clocs(end);
    
    ClocsNTend = RF-1 + 3*Clocs;
    
    ConsSeq(BctCons95==0) = 'X';
    ConservedSeq{y} = ConsSeq(1:ClocsNTend);
end

%Now, find the C locs of all Vmap
for j = 1:size(Vmap,1)
    VmapFamNum = Vmap{j,4};
    UnqFamNumLoc = findHeader(UnqFamNum,VmapFamNum);
    RefSeq = ConservedSeq{UnqFamNumLoc}(end-30:end);
    CurSeq = Vmap{j,1};
    [~, Alignment, StartAt, ~] = convolveSeq(CurSeq,RefSeq,0,0);
    Cloc = abs(StartAt(2))+length(RefSeq)-2; %Location from left to right
    ClocDist = length(CurSeq) - Cloc + 1; %Distance from end until left side of TGT seq. (end-ClocDist+1:end-ClocDist+3) should give you TGT/.   
    
    if abs(ClocDist) > 30 || ClocDist < 0 %Hard to get TGT greater than 30
        ClocDist = 0;
        disp('Pseudogene, undefined C');
    end
    Vmap{j,end} = ClocDist;
end
