%findCDR will fill in the CDR1, CDR2, and CDR3 information, such as AA, AA
%length, and start and end DNA location of this CDR for each sequence entry
%in VDJdata.
%
%  VDJdata = findCDR(VDJdata, Map, DB, CDRNum)
%
%  VDJdata = findCDR(VDJdata, Map, DB, ..., 'imgt')
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    Map: structure of indeix BRILIA data columns
%    DB: Gene database structure (getCurrentDatabase.m)
%    CDRNum: Vector of CDR number(s). Default is [1:3]
%    'imgt': Use IMGT's CDR3 position, which EXCLUDES the 104 & 118 AA.
%            Otherwise, BRILIA uses the 104/118 as anchors for annotation.
%            NOTE: this only is used for CDR3 (CDRNum == 3);
%
%  OUTPUT
%    VDJdata: modified VDJdata with the CDR information filled
%
%  NOTE
%    This must be run after finding V(D)J annotations.
%
function VDJdata = findCDR(VDJdata, Map, DB, varargin)
if isempty(VDJdata); return; end
UseIMGTLoc = strcmpi(varargin, 'IMGT');
UseIMGT = any(UseIMGTLoc);
varargin = varargin(~UseIMGTLoc);

if isempty(varargin)
    Num = 1:3;
else
    Num = varargin{1};
end

for n = 1:length(Num)
    for c = 1:length(Map.Chain)
        C = lower(Map.Chain(c));
        switch Num(n)
            case 1
                VDJdata = findCDR_Prefilter(VDJdata, Map, DB, C, Num(n));
            case 2
                VDJdata = findCDR_Prefilter(VDJdata, Map, DB, C, Num(n));
            case 3
                VDJdata = findCDR3_Prefilter(VDJdata, Map, DB, C, Num(n), UseIMGT);
            otherwise
                error('%s: Num must be a vector between 1:3', mfilename);
        end
    end
end

function VDJdata = findCDR_Prefilter(VDJdata, Map, DB, C, Num)
DelIdx = Map.([C 'Del'])(1);
SeqIdx = Map.([C 'Seq']);
LenIdx = Map.([C 'Length'])(1);
GeneIdx = Map.([C 'GeneName'])(1);KeepLoc = ~any(cellfun('isempty', VDJdata(:, [DelIdx; SeqIdx; LenIdx; GeneIdx])), 2);

if strcmpi(C, 'H')
    VDJdata(KeepLoc, :) = findCDR_Calc(VDJdata(KeepLoc, :), Map, DB, C, 'V', Num);
else 
    KLoc = KeepLoc & contains(VDJdata(:, GeneIdx), 'IGKV', 'ignorecase', true);
    VDJdata(KLoc, :) = findCDR_Calc(VDJdata(KLoc, :), Map, DB, C, 'Vk', Num);
    
    LLoc = KeepLoc & contains(VDJdata(:, GeneIdx), 'IGLV', 'ignorecase', true);
    VDJdata(LLoc, :) = findCDR_Calc(VDJdata(LLoc, :), Map, DB, C, 'Vl', Num);
end

function VDJdata = findCDR_Calc(VDJdata, Map, DB, C, X, Num)
DelIdx = Map.([C 'Del'])(1);
SeqIdx = Map.([C 'Seq']);
LenIdx = Map.([C 'Length'])(1);
GeneIdx = Map.([C 'GeneName'])(1);
CDRIdx = Map.([C 'CDR' num2str(Num)]);
DBIdx = getGeneNameIdx(VDJdata(:, GeneIdx), DB, X);
UnqDBIdx = unique(DBIdx);
M = getMapHeaderVar(DB.MapHeader);
for j = 1:length(UnqDBIdx)
    if UnqDBIdx(j) == 0; continue; end
    RefVLen = length(DB.([X 'map']){UnqDBIdx(j), M.Seq});
    RefCDRS = DB.([X 'map']){UnqDBIdx(j), M.(sprintf('CDR%ds', Num))};
    RefCDRE = DB.([X 'map']){UnqDBIdx(j), M.(sprintf('CDR%de', Num))};
    RefCDR  = nt2aa(DB.([X 'map']){UnqDBIdx(j), M.Seq}(RefCDRS:RefCDRE), 'ACGTonly', false, 'alternativestart', false);

    Idx = find(UnqDBIdx(j) == DBIdx);
    SeqVLen = cell2mat(VDJdata(Idx, LenIdx));
    SeqVDel = cell2mat(VDJdata(Idx, DelIdx));
    
    Vshift = SeqVLen + SeqVDel - RefVLen;
    CDRS = Vshift + RefCDRS;
    CDRE = Vshift + RefCDRE;
    CDRLen = (CDRE - CDRS + 1) / 3;
    
    SeqCDR = cell(size(Idx));
    for k = 1:length(Idx)
        if min([CDRS(k) CDRE(k)]) > 0
            SeqCDR{k} = nt2aa(VDJdata{Idx(k), SeqIdx}(CDRS(k):CDRE(k)), 'ACGTonly', false, 'alternativestart', false);
            MissLoc = SeqCDR{k} ~= RefCDR;
            SeqCDR{k}(MissLoc) = lower(SeqCDR{k}(MissLoc));
        else
            SeqCDR{k} = RefCDR;
        end
    end
    VDJdata(Idx, CDRIdx) = [SeqCDR num2cell([CDRLen CDRS CDRE])];
end
        
function VDJdata = findCDR3_Prefilter(VDJdata, Map, DB, C, Num, UseIMGT)
DelIdx = Map.([C 'Del'])([1 end]);
SeqIdx = Map.([C 'Seq']);
LenIdx = Map.([C 'Length'])([1 end]);
GeneIdx = Map.([C 'GeneName'])([1 end]);
KeepLoc = ~any(cellfun('isempty', VDJdata(:, [DelIdx; SeqIdx; LenIdx; GeneIdx])), 2);

if Num < 3
    if strcmpi(C, 'H')
        VDJdata(KeepLoc, :) = findCDR_Calc(VDJdata(KeepLoc, :), Map, DB, C, 'V', Num);
    else 
        KLoc = KeepLoc & contains(VDJdata(:, GeneIdx), 'IGKV', 'ignorecase', true);
        VDJdata(KLoc, :) = findCDR_Calc(VDJdata(KLoc, :), Map, DB, C, 'Vk', Num);

        LLoc = KeepLoc & contains(VDJdata(:, GeneIdx), 'IGLV', 'ignorecase', true);
        VDJdata(LLoc, :) = findCDR_Calc(VDJdata(LLoc, :), Map, DB, C, 'Vl', Num);
    end
else
    if strcmpi(C, 'H')
        VDJdata(KeepLoc, :) = findCDR3_Calc(VDJdata(KeepLoc, :), Map, DB, C, 'V', UseIMGT);
    else 
        KLoc = KeepLoc & contains(VDJdata(:, GeneIdx), 'IGKV', 'ignorecase', true);
        VDJdata(KLoc, :) = findCDR3_Calc(VDJdata(KLoc, :), Map, DB, C, 'Vk', UseIMGT);

        LLoc = KeepLoc & contains(VDJdata(:, GeneIdx), 'IGLV', 'ignorecase', true);
        VDJdata(LLoc, :) = findCDR3_Calc(VDJdata(LLoc, :), Map, DB, C, 'Vl', UseIMGT);
    end
end

function VDJdata = findCDR3_Calc(VDJdata, Map, DB, C, V, UseIMGT)
if isempty(VDJdata); return; end
DelIdx = Map.([C 'Del'])([1 end]);
SeqIdx = Map.([C 'Seq']);
LenIdx = Map.([C 'Length'])([1 end]);
GeneIdx = Map.([C 'GeneName'])([1 end]);
CDRIdx = Map.([C 'CDR3']);
switch V
    case 'V'
        J = 'J';
    case 'Vk'
        J = 'Jk';
    case 'Vl'
        J = 'Jl';
end
VDBIdx = getGeneNameIdx(VDJdata(:, GeneIdx(1)),   DB, V);
JDBIdx = getGeneNameIdx(VDJdata(:, GeneIdx(end)), DB, J);

UnqVDBIdx = unique(VDBIdx);
M = getMapHeaderVar(DB.MapHeader);
for j = 1:length(UnqVDBIdx)
    if UnqVDBIdx(j) == 0; continue; end
    Idx = find(UnqVDBIdx(j) == VDBIdx);
    
    RefVAnchor = DB.([V 'map']){UnqVDBIdx(j), M.Anchor};
    RefJAnchor = vertcat(DB.([J 'map']){JDBIdx(Idx), M.Anchor});

    SeqLen = cellfun('length', VDJdata(Idx, SeqIdx));
    SeqVLen = cell2mat(VDJdata(Idx, LenIdx(1)));
    SeqVDel = cell2mat(VDJdata(Idx, DelIdx(1)));
    SeqJLen = cell2mat(VDJdata(Idx, LenIdx(end)));
    SeqJDel = cell2mat(VDJdata(Idx, DelIdx(end)));
    
    CDRS = SeqVLen + SeqVDel - RefVAnchor + 1; %Add 1 to ensure it starts on that nt
    CDRE = SeqLen - (SeqJLen + SeqJDel) + RefJAnchor + 2; %Add 2 to ensure it end on 3rd nt of codon   
    if UseIMGT
        CDRS = CDRS + 3;
        CDRE = CDRE - 3;
    end
    CDRLen = (CDRE - CDRS + 1) / 3;
    
    SeqCDR = cell(size(Idx));
    for k = 1:length(Idx)
        if min([CDRS(k) CDRE(k)]) > 0 && CDRE(k) <= length(VDJdata{Idx(k), SeqIdx})
            SeqCDR{k} = nt2aa(VDJdata{Idx(k), SeqIdx}(CDRS(k):CDRE(k)), 'ACGTonly', false, 'alternativestart', false);
        end
    end
    
    VDJdata(Idx, CDRIdx) = [SeqCDR num2cell([CDRLen CDRS CDRE])];
end