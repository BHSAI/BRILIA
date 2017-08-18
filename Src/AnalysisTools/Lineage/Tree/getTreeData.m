%getTreeData2 will extract necessary information required to plot a tree
%given VDJdata from ONE cluster (same group number).
%
%  INPUT
%    Tdata: the MxN data subset of VDJdata for the same group number. Will
%      complain if multiple group numbers are detected.
%    VDJheader: the 1xN cell of data header names of VDJdata
%
%  OUPUT
%    AncMapStruct: Structure holding various Mx3 AncMap
%      .HAM = AncMap with par2child hamming distance
%      .SHM = AncMap with par2child BRILIA shm distance
%      .HAMPERC = AncMap with par2child hamming distance / seq length
%      .SHMPERC = AncMap with par2child BRILIA shm distance / seq length
%    TreeName: suggested name of the tree based on VDJ and VJ annotation.
%      If H chain, {'Vgene | Dgene | Jgene' 'Grp #, Size #, TC #'}
%      If L chain, {'Vxgene | Jxgene' 'Grp #, Size #, TC #'}
%      If HL chains, {'Vgene | Dgene | Jgene' 'Vxgene | Jxgene' 
%                                                   'Grp #, Size #, TC #'}
%      where Size of the group in number of unique sequence, TC is total
%      template count of the group.
%    CDR3Names: is a Mx1 cell of CDR3 names used for making color legends.
%      IF H or L chain, CDR3Names is just CDR3's
%      If HL chain, CDR3Names is joined by a dash as CDR3H-CDR3L
%    TemplateCount: Template count to use for drawing tree. For germline, 
%      will set to 0.
%
%  NOTE
%    IF germline is the same as the 1st sequence of the group, will return
%    all values of the same number of entries as Tdata.
%    IF germline is different from the 1st sequence, will add the germline
%    sequence as an entry but with 0 template count, returning parameters
%    with size(Tdata,1) + 1 entries.

function [AncMapS, TreeName, CDR3Name, TemplateCount] = getTreeData(Tdata, VDJheader)
AncMapS = [];
TreeName = [];
CDR3Name = [];

%Check if there are multiple groups
[H, L, Chain] = getAllHeaderVar(VDJheader);
GrpNum = cell2mat(Tdata(:, H.GrpNumLoc));
UnqGrpNum = unique(GrpNum);
if length(UnqGrpNum) > 1
    warning('%s: Tdata must contain annotations from the same cluster. No tree was made.', mfilename);
    return;
end

%Determine key tree parameter data
SeqNumLoc = H.SeqNumLoc;
TemplateLoc = H.TemplateLoc;
CDR3Locs = [];
CDR3sLocs = [];
CDR3eLocs = [];
SeqLocs = [];
RefSeqLocs = [];
GeneNameLocs = {};
for j = 1:length(Chain)
    if strcmpi(Chain(j), 'H')
        B = H;
    else
        B = L;
    end
    try
        CDR3Locs = cat(2, CDR3Locs, B.CDR3Loc(1));
        CDR3sLocs = cat(2, CDR3sLocs, B.CDR3Loc(3));
        CDR3eLocs = cat(2, CDR3eLocs, B.CDR3Loc(4));
        SeqLocs = cat(2, SeqLocs, B.SeqLoc);
        RefSeqLocs = cat(2, RefSeqLocs, B.RefSeqLoc);
        GeneNameLocs = cat(2, GeneNameLocs, B.GeneNameLoc);
    catch
    end
end

%Make sure all values are present
CDR3Locs(CDR3Locs == 0) = [];
CDR3sLocs(CDR3sLocs == 0) = [];
CDR3eLocs(CDR3eLocs == 0) = [];
SeqLocs(SeqLocs == 0) = [];
RefSeqLocs(RefSeqLocs == 0) = [];

%--------------------------------------------------------------------------
%Getting TreeName

TreeName = cell(length(GeneNameLocs) + 1, 1);
for j = 1:length(GeneNameLocs)
    GeneNames = Tdata(1, GeneNameLocs{j});
    for w = 1:length(GeneNames)
        TempNames = regexpi(GeneNames{w},'\|','Split'); %Get the first recommended name
        TempName = strrep(TempNames{1}, 'IGH', ''); %Remove IGH for heavy chain
        GeneNames{w} = strrep(TempName, 'IG', ''); %Remove IG for light chain. Save.
    end
    RepPat = repmat(' %s |',1,length(GeneNames));
    RepPat([1 end]) = [];
    TreeName{j} = sprintf(RepPat, GeneNames{:});
end
TreeName{end} = sprintf('Grp %d, Size %d, TC %d', UnqGrpNum, size(Tdata, 1), sum(cell2mat(Tdata(:, TemplateLoc))));

%--------------------------------------------------------------------------
%Getting AncMapS

AncMapAll = zeros(size(Tdata, 1), 7); %[Child Par HAMdist SHMdist HAMperc SHMperc TemplateCount];
for j = 1:size(Tdata, 1)
    %Determine this child seq's absolute parent seq number
    ChildSeqNum = Tdata{j, SeqNumLoc};
    ParentSeqNum = 0;
    ParentSeq = Tdata(j, RefSeqLocs);
    for g = 1:size(Tdata, 1)
        if g == j; continue; end
        Match = zeros(1, length(SeqLocs));
        for k = 1:length(SeqLocs)
            if strcmpi(ParentSeq{k}, Tdata{g, SeqLocs(k)})
                Match(k) = 1;
            end
        end
        if min(Match) == 1
            ParentSeqNum = Tdata{g, SeqNumLoc};
            break;
        end
    end

    %Determine the SHMdist and HAMdist and template count
    SHMdist = 0; %Parent to child SHM distance
    HAMdist = 0; %Parent to child HAM distance
    SeqLen = 0;
    for k = 1:length(SeqLocs)
        ChildSeq = Tdata{j, SeqLocs(k)};
        ParentSeq = Tdata{j, RefSeqLocs(k)};
        [SHMdistT, ~, HAMdistT] = calcSHMHAMdist(ParentSeq, ChildSeq);
        SHMdist = SHMdist + SHMdistT;
        HAMdist = HAMdist + HAMdistT;
        SeqLen = SeqLen + length(ParentSeq);
    end
    SHMperc = SHMdist / SeqLen * 100; 
    HAMperc = HAMdist / SeqLen * 100;

    %Determine Template Count
    TemplateCount = Tdata{j, TemplateLoc};
    if isempty(TemplateCount) || TemplateCount <= 0
        TemplateCount = 1;
    end
    
    AncMapAll(j, :) = [ChildSeqNum ParentSeqNum HAMdist SHMdist HAMperc SHMperc TemplateCount];
end

%Determine if you need to add the germline with template count 0
AddedGermline = 0; %To tell CDR3Name code to fetch germline CDR3
if AncMapAll(1, 2) == 0 && AncMapAll(1, 3) > 0
    AncMapAll = cat(1, zeros(1, size(AncMapAll, 2)), AncMapAll);
    AncMapAll(:, 1:2) = AncMapAll(:, 1:2) + 1;
    AncMapAll(2, 2) = 1; %Link 1st (now 2nd) to germline (now 1st)
    AncMapAll(1, 2) = 0; %Ensure it's always 0.
    AddedGermline = 1;
end
AncMapAll = renumberAncMap(AncMapAll);

%Create the structured version of AncMap
DistOrder = {'HAM', 'SHM', 'HAMPERC', 'SHMPERC'};
for k = 1:length(DistOrder)
    AncMapS.(DistOrder{k}) = AncMapAll(:, [1:2 k+2 size(AncMapAll, 2)]);
end

%--------------------------------------------------------------------------
%Getting CDR3Name

CDR3Name = cell(size(Tdata, 1) + AddedGermline, 1);
for j = 1:size(Tdata, 1)
    %Combine and extract CDR3 Info
    if length(CDR3Locs) > 1
        RepPat = repmat('%s:', 1, length(CDR3Locs));        
        CDR3seq = sprintf(RepPat, Tdata{j, CDR3Locs});
        CDR3seq(end) = [];
        CDR3Name{j + AddedGermline} = CDR3seq;
    else
        CDR3Name{j + AddedGermline} = Tdata{j, CDR3Locs};
    end
end
    
if AddedGermline == 1
    %Add in the germline sequence CDR3
    CDR3NameGerm = cell(1, length(CDR3Locs));
    for k = 1:length(CDR3Locs)
        if ~isempty(CDR3sLocs) && ~isempty(CDR3eLocs)
            RefSeq = Tdata{1, RefSeqLocs(k)};
            CDR3s = Tdata{1, CDR3sLocs(k)};
            CDR3e = Tdata{1, CDR3eLocs(k)};
            CDR3NameGerm{k} = nt2aa(RefSeq(CDR3s:CDR3e), 'ACGTonly', 'false');
        end
    end
    
    %Combine and extract CDR3 Info
    RepPat = repmat('%s:', 1, length(CDR3NameGerm));        
    CDR3seq = sprintf(RepPat, CDR3NameGerm{:});
    CDR3seq(end) = [];
    CDR3Name{1} = CDR3seq;
end

%--------------------------------------------------------------------------
%Getting TemplateCount

TemplateCount = zeros(size(Tdata, 1) + AddedGermline, 1);
TemplateCount(1+AddedGermline:end) = cell2mat(Tdata(:, TemplateLoc));
