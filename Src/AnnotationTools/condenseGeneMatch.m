%condenseGeneMatch is used by findGeneMatch to reduce VDJ matches with the same
%gene family, gene deletion, and alignment score. This is used to "reduce"
%the number of degenerate match outputs by condensing the family names into
%something like: "IGHV01: IGHV01-01*01, IGHV01-01*02". Will return all
%unique gene families. Selecting which one of the degenerate solution is
%best must be determined by a later function or by findGeneMatch.
%
%  ReducedMatch = condenseGeneMatch(GeneMatch)
%
%  INPUT
%    GeneMatch: Mx6 cell matrix output from findGeneMatch
%      Col1   Gene number in Xmap
%      Col2   Full gene name(s)
%      Col3   [LeftLength MiddleLength RighLength] of the reference gene
%      Col4   [LeftLength MiddleLength RighLength] of the Seq
%      Col5   [(# of matches) AlignmentScore]
%      Col6   3xN character alignment results
%
%  OUTPUT
%    ReducedMatch: Mx6 cell matrix output like from findGeneMatch, but
%      gene matches with the same alignment score, and Left, Middle, and
%      Right lengths of the gene seq and sample seq are grouped together.
%
%  EXAMPLE
%    GeneMatch = {1 'IGHV1-01' [40 100 3] [0 100 20] [90 450] 'ACGTG';
%                 2 'IGHV1-02' [40 100 3] [0 100 20] [90 450] 'ACGTG';
%                 3 'IGHV1-03' [41 101 2] [0  90 30] [30 350] 'AGGTG'};
%    ReducedMatch = condenseGeneMatch(GeneMatch)
%    ReducedMatch = 
%       [3]    'IGHV1-03'          [1x3 doub] [1x3 doub] [1x2 doub] 'AGGTG'
%       [1 2]  'IGHV1-01|IGHV1-02' [1x3 doub] [1x3 doub] [1x2 doub] 'ACGTG'
%
%  See also findGeneMatch, findVDJmatch, find VJmatch

function ReducedMatch = condenseGeneMatch(GeneMatch)
%Determine entries with same annotations and alignment scores
RefLMR = cell2mat(GeneMatch(:,3)); %Ref gene Left Mid Right lengths
SamLMR = cell2mat(GeneMatch(:,4)); %Sample gene Left Mid Right lengths
Scores = cell2mat(GeneMatch(:,5)); %Alignment scores

%Find those with the same scores and LMRs, which will be grouped
SortMat = [Scores(:,end) SamLMR RefLMR]; %Use the conse alignment score, if available
[~,~,UnqIdx] = unique(SortMat,'rows');

%Condense GeneMatch per unique SortMat row entity.
ReducedMatch = cell(max(UnqIdx),size(GeneMatch,2));
j = 1;
for y = 1:max(UnqIdx)    
    %Condense the gene names
    Idx = find(UnqIdx == y);
    RepPat = repmat('%s|',1,length(Idx));
    RepPat(end) = [];
    NewName = sprintf(RepPat,GeneMatch{Idx,2});
    
    %Update the UnqGeneMatch
    ReducedMatch(j,:) = GeneMatch(Idx(1),:);
    ReducedMatch{j,1} = cell2mat(GeneMatch(Idx,1))';
    ReducedMatch{j,2} = NewName;
    
    j = j+1;
end   
