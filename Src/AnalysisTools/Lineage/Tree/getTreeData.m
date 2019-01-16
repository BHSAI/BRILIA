%getTreeData will extract necessary information required to plot a tree
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

function [AncMapS, TreeName, CDR3Name, TemplateCount] = getTreeData(Tdata, Map)
AncMapS = [];
TreeName = [];
CDR3Name = [];
TemplateCount = 0;

Map = getVDJmapper(Map);

%Check if there are multiple groups
GrpNum = cell2mat(Tdata(:, Map.GrpNum));
UnqGrpNum = unique(GrpNum);
if length(UnqGrpNum) > 1
    warning('%s: Tdata must contain annotations from the same cluster. No tree was made.', mfilename);
    return
end

%Determine key tree parameter data
SeqNumIdx = Map.SeqNum;
ParNumIdx = Map.ParNum;
TemplateIdx = Map.Template;
SeqIdx = nonzeros([Map.hSeq; Map.lSeq]);
RefSeqIdx = nonzeros([Map.hRefSeq; Map.lRefSeq]);
CDR3Idx = zeros(length(Map.Chain), 1);
CDR3SIdx = zeros(length(Map.Chain), 1);
CDR3EIdx = zeros(length(Map.Chain), 1);
GeneNameIdx = cell(length(Map.Chain), 1);
for c = 1:length(Map.Chain)
    C = lower(Map.Chain(c));
    CDR3Idx(c) = Map.([C 'CDR3'])(1);
    CDR3SIdx(c) = Map.([C 'CDR3'])(3);
    CDR3EIdx(c) = Map.([C 'CDR3'])(4);
    GeneNameIdx{c} = Map.([C 'GeneName']);
end

AncMap = zeros(size(Tdata, 1), 5); %[Child Par HAMdist HAMperc TemplateCount];
AncMap(:, [1:2 5]) = cell2mat(Tdata(:, [SeqNumIdx; ParNumIdx; TemplateIdx]));
AncMap = renumberAncMap(AncMap);
if any(findTreeCycle(AncMap))
    error('%s: Found a cyclic dependency in group number %d.', mfilename, UnqGrpNum);
end

for j = 1:size(Tdata, 1)
    MatchCt = 0;
    SeqLen = 0;
    if AncMap(j, 2) == 0 %RefSeq is the inferred germline
        RPos = j; 
        RIdx = RefSeqIdx;
    else %RefSeq is another sequence
        RPos = AncMap(j, 2);
        RIdx = SeqIdx;
    end
    for c = 1:length(SeqIdx)
        MatchCt = MatchCt + sum(cmprSeqMEX(Tdata{j, SeqIdx(c)}, Tdata{RPos, RIdx(c)}, 'n'));
        SeqLen = SeqLen + length(Tdata{RPos, RIdx(c)});
    end
    AncMap(j, 3:4) = [(SeqLen-MatchCt) ((SeqLen-MatchCt)/SeqLen*100)]; %HamDist, HamPerc
end

%Add a germline sequence in case the 1st observed sequence has HamDist > 0 to germline
NeedGermline = AncMap(1, 2) == 0 && AncMap(1, 3) > 0; %To tell CDR3Name code to fetch germline CDR3
if NeedGermline
    AncMap = cat(1, zeros(1, size(AncMap, 2)), AncMap);
    AncMap(:, 1:2) = AncMap(:, 1:2) + 1;
    AncMap(2, 2) = 1; %Link 1st (now 2nd) to germline (now 1st)
    AncMap(1, 2) = 0; %Ensure it's always 0.
end
AncMap = renumberAncMap(AncMap);

%Create the structured version of AncMap
DistOrder = {'HAM', 'HAMPERC'};
for k = 1:length(DistOrder)
    AncMapS.(DistOrder{k}) = AncMap(:, [1:2 k+2 size(AncMap, 2)]);
end

%--------------------------------------------------------------------------
%Getting TreeName
if nargout >= 2
    TreeName = cell(length(GeneNameIdx) + 1, 1);
    for j = 1:length(GeneNameIdx)
        GeneNames = Tdata(1, GeneNameIdx{j});
        for w = 1:length(GeneNames)
            TempNames = strsplit(GeneNames{w}, '|');  %Get the first recommended name
            GeneNames{w} = strrep(TempNames{1}, 'IG', ''); 
        end
        RepPat = repmat(' %s |',1,length(GeneNames));
        RepPat([1 end]) = [];
        TreeName{j} = sprintf(RepPat, GeneNames{:});
    end
    TreeName{end} = sprintf('Grp %d, Size %d, TC %d', UnqGrpNum, size(Tdata, 1), sum(cell2mat(Tdata(:, TemplateIdx))));
end

%--------------------------------------------------------------------------
%Getting CDR3Name
if nargout > 3
    CDR3Name = cell(size(Tdata, 1) + NeedGermline, 1);
    for j = 1:size(Tdata, 1)
        %Combine and extract CDR3 Info
        if length(CDR3Idx) > 1
            RepPat = repmat('%s:', 1, length(CDR3Idx));        
            CDR3seq = sprintf(RepPat, Tdata{j, CDR3Idx});
            CDR3seq(end) = [];
            CDR3Name{j + NeedGermline} = CDR3seq;
        else
            CDR3Name{j + NeedGermline} = Tdata{j, CDR3Idx};
        end
    end

    if NeedGermline
        %Add in the germline sequence CDR3
        CDR3NameGerm = cell(1, length(CDR3Idx));
        for k = 1:length(CDR3Idx)
            if ~isempty(CDR3SIdx) && ~isempty(CDR3EIdx)
                RefSeq = Tdata{1, RefSeqIdx(k)};
                CDR3s = Tdata{1, CDR3SIdx(k)};
                CDR3e = Tdata{1, CDR3EIdx(k)};
                CDR3NameGerm{k} = nt2aa(RefSeq(CDR3s:CDR3e), 'ACGTonly', false);        
            end
        end
        

        %Combine and extract CDR3 Info
        CDR3seq = sprintf('%s:', CDR3NameGerm{:});
        CDR3Name{1} = CDR3seq(1:end-1);
        
        if contains(CDR3Name{1}, '*')
            CDR3Name{1} = CDR3Name{2}; %Taking ancestral sequence instead
        end
    end
end

%--------------------------------------------------------------------------
%Getting TemplateCount after all the germline edits
if nargout >= 4
    TemplateCount = AncMap(:, 5);
end 
