%trimMatchResults will take a binary array of sequence matches and trim
%matches (1's) from the left, right, or both direction(s) until it
%encounters a sequence block where 3 out of 4 are matches. Used mainly for
%convolveSeq.
%
%  MatchResults = trimMatchRestuls(MatchResults,TrimSide) 
%
%  [MatchResults,TrimIdx] = trimMatchRestuls(MatchResults,TrimSide)
%
%  INPUT
%    MatchResults: MxN logical matrix of sequence match. 
%    TrimSide ['left','right','both']: Trims this side of the MatchResults.
%
%  OUPUT
%    MatchResults: MxN logical matrix of trimmed MatchResults.
%    TrimIdx: MxN logical matrix showing trimmed sites.
%
%  NOTES 
%    If you are allowing misses, such as [1 0 1], to be counted as matches,
%    then fill in those allowed miss with 1's before this.
%
%  EXAMPLE
%    MatchResults(1,:) = [1 1 0 0 1 0 0 0 1 1 1 0 0 1]>0;
%    MatchResults(2,:) = [1 1 1 0 0 1 0 1 0 1 0 1 0 1]>0;
%
%    TrimmedLeft = trimMatchResults(MatchResults,'left')
%    TrimmedLeft =
%      0    0    0    0    0    0    0    0    1    1    1    0    0    1
%      1    1    1    0    0    1    0    1    0    1    0    1    0    1
%      
%    TrimmedRight = trimMatchResults(MatchResults,'right')
%    TrimmedRight =
%      1    1    0    0    1    0    0    0    1    1    1    0    0    0
%      1    1    1    0    0    0    0    0    0    0    0    0    0    0
%
%    TrimmedBoth = trimMatchResults(MatchResults,'both')
%    TrimmedBoth =
%      0    0    0    0    0    0    0    0    1    1    1    0    0    0
%      1    1    1    0    0    0    0    0    0    0    0    0    0    0
%
%  See also makeDiagonalSeq, convolveSeq

function [MatchResults,varargout] = trimMatchResults(MatchResults,TrimSide)
%Input checking and formatting
if ~islogical(MatchResults)
    error('Error in %s: MatchResults must be a logical matrix.',mfilename);
end
if ~ismember(lower(TrimSide),{'both','left','right'})
    error('Error in %s: TrimSide not defined correctly.',mfilename);
end
TrimSide = lower(TrimSide(1));

%Tracks where trims were done. Needed for calcAlignScore.
TrimIdx = zeros(size(MatchResults),'logical'); 

%If the matrix is too small, just skip this entirely.
if size(MatchResults,2) <= 4
    if nargout == 2;
        varargout{1} = TrimIdx;
    end
    return;
end
     
%Determine where there are 3/4 matches are, marked as ValidLoc.
SumMat = MatchResults(:,1:end-3);
for j = 2:4
    SumMat = SumMat + MatchResults(:,j:end-4+j);
end
ValidLoc = SumMat >= 3;

%Modify MatchResults and TrimIdx to reflect trimmed locations.
for j = 1:size(SumMat)
    ValidIdx = find(ValidLoc(j,:));
    
    %No match case
    if isempty(ValidIdx);
        TrimIdx(j,:) = 1;
        MatchResults(j,:) = 0;
        continue
    end
    
    %Do left side trimming
    if TrimSide == 'b' || TrimSide == 'l'
        MatchResults(j,1:ValidIdx(1)-1) = 0;
        TrimIdx(j,1:ValidIdx(1)) = 1;
        TrimIdx(j,ValidIdx(1)) = ~MatchResults(j,ValidIdx(1));
    end
    
    %Do right side trimming
    if TrimSide == 'b' || TrimSide == 'r'
        MatchResults(j,ValidIdx(end)+4:end) = 0;
        TrimIdx(j,ValidIdx(end)+3:end) = 1;
        TrimIdx(j,ValidIdx(end)+3) = ~MatchResults(j,ValidIdx(end)+3);
    end
end

%Return TrimIdx if summoned
if nargout == 2
    varargout{1} = TrimIdx;
end
