%getTreeData will go through VDJdata and extract the parent-child
%relationship table. This assumes parent-child relationships are stored
%properly and created by clusterGene.m. This code produces TreeData that is
%minimally needed to plot lineage trees.
%
%  [TreeData, TreeHeader] = getTreeData
%
%  [TreeData, TreeHeader] = getTreeData(FileName)
%
%  [TreeData, TreeHeader] = getTreeData(VDJdata,NewHeader)
%
%  INPUT
%    VDJdata: MxN cell of VDJ annotations AFTER clustering by lineage
%    NewHeader: 1xN cell storing header names for VDJdata
%
%  OUTPUT
%    TreeData: a Mx5 matrix required to draw a lineage tree. The column
%    information are as follows:
%      1) Current sequence number
%      2) Parent sequence number
%      3) Template count of current sequence
%      4) SHM distance from parent to current sequence
%      5) Hamming distance from parent to current sequence
%      6) CDR3 sequence
%      7) Group name as Vname|Dname|Jname (Grp##, Size = ##)
%      8) Group number
%
%  NOTE
%    If no template count data are provided, will assume each sequence has
%    a template count of 1.
%
%    The 1st sequence in every new group is assumed most ancestral and is
%    linked to the germline sequence. Hence the parent sequence number is
%    0.
%
%    To plot a lineage tree, use plotTreeData.
%
%    To determine lineage relationships, use buildTreeData.
%  
%  See also plotTreeData, buildTreeData

function [TreeData, TreeHeader] = getTreeData(varargin)
%See if user gave VDJdata
if isempty(varargin) || (~isempty(varargin) && isempty(varargin{1})) %Need to find file
    [VDJdata,NewHeader] = openSeqData;
elseif ~isempty(varargin)
    if ischar(varargin{1}) %Filename was given
        [VDJdata,NewHeader] = openSeqData(varargin{1});
    elseif length(varargin) == 2 && iscell(varargin{1}) && iscell(varargin{2}) %VDJdata and NewHeader was given
        VDJdata = varargin{1};
        NewHeader = varargin{2};
    end
else
    error('getTreeData: Check the inputs');
end
getHeaderVar;

%Extract sequence grouping information
GrpNum = cell2mat(VDJdata(:,GrpNumLoc));
UnqGrpNum = unique(GrpNum);

%Initialize output variables
TreeHeader = {'ChildSeqNum' 'ParentSeqNum' 'TemplateCount' 'SHMdist' 'HAMdist' 'CDR3seq' 'GrpName' 'GrpNum'};
TreeData = cell(size(VDJdata,1),length(TreeHeader)); 
TreeData(:,1) = num2cell(GrpNum); %Save the group numbers

%Fill in TreeData group by group
for y = 1:length(UnqGrpNum)
    GrpIdx = find(UnqGrpNum(y) == GrpNum);

    %Build the default group name as Vname|Dname|Jname (Grp##, Size = ##)
    VDJnames = VDJdata(GrpIdx(1),FamLoc);
    for w = 1:3
        TempName = VDJnames{w};
        ColonLoc = regexp(TempName,'\:');
        if isempty(ColonLoc)
            ColonLoc = 0;
        end
        TempName = TempName(ColonLoc+1:end);
        TempName = strrep(TempName,' ','');            
        TempName = regexp(TempName,'\|','split');
        TempName = strrep(TempName{1},'IGH','');

        %Remove the "0" in V01, D01, etc.
        if TempName(2) == '0'
            TempName(2) = [];
        end
        VDJnames{w} = TempName;
    end
    GrpName = sprintf('%s | %s | %s (Grp %d, Size %d)',VDJnames{1},VDJnames{2},VDJnames{3},UnqGrpNum(y),length(GrpIdx));
        
    for j = 1:length(GrpIdx)
        %The current Seq is the child sequence
        ChildSeqNum = VDJdata{GrpIdx(j),SeqNumLoc};
        ChildRefSeq = VDJdata{GrpIdx(j),RefSeqLoc}; %The parent seq matches this child's RefSeq.
        
        %The parent of current Seq is the RefSeq, which must be found.
        if j == 1 %1st sequence is always linked to germline sequence
            ParentSeqNum = 0;
        else
            GrpName = ''; %Don't need it anymore
            
            %Find the nth GrpIdx position that has Seq = ChildRefSeq.
            ParentSeqRelIdx = 0;
            for k = 1:length(GrpIdx)
                if strcmp(VDJdata{GrpIdx(k),SeqLoc},ChildRefSeq) == 1
                    ParentSeqRelIdx = k;
                    break
                end
            end
            
            if ParentSeqRelIdx == 0 %No parent found, defaulting to germline.
                warning('Cannot resolve parent of sequence #%d. Defaulting to germline seq.',ChildSeqNum);
                ParentSeqNum =  0;
            else %Return the sequence number
                ParentSeqNum = VDJdata{GrpIdx(ParentSeqRelIdx(1)),SeqNumLoc};
            end
        end
        
        %Determine Template Count        
        TemplateCount = VDJdata{j,TemplateLoc};
        if isempty(TemplateCount) || TemplateCount <= 0
            TemplateCount = 1;
        end
        
        %Determine the SHMdist and HAMdist and template count
        ChildSeq = VDJdata{GrpIdx(j),SeqLoc};
        ParentSeq = VDJdata{GrpIdx(j),RefSeqLoc};
        [SHMdist, ~, HAMdist] = calcSHMHAMdist(ParentSeq,ChildSeq);
        
        %Extract CDR3info, if available
        CDR3seq = VDJdata{GrpIdx(j),CDR3Loc(1)};
        
        %Fill in TreeData
        TreeData(GrpIdx(j),:) = {ChildSeqNum ParentSeqNum TemplateCount SHMdist HAMdist CDR3seq GrpName UnqGrpNum(y)};
    end
end