%flipNregion will analyze a sequence in the N region, and determine if it
%was made from TDT 3' end of coding or noncoding. Returns the complement
%strand IF there are more nucleotides that are assocaited with the
%complement synthesis.
%
%Note: The assumption is that 3' additions will add G's and A's. If it's on
%the leading strand, you'll see G's and A's, where as if it's on the other
%strand, you'll get C's and T's. 
%
%  Examples:
%  Vague scenarios (no flip can be decided. Need more information.)
%    LeftNTs = 'AG'
%    RightNTs = 'CT'
%    Seq1 = 'AGGGGCT'
%    Seq2 = 'AGCCCCT'
%    [Seq1TDT, Flipped] = flipNregion(Seq1,LeftNTs,RightNTs)
%      Seq1TDT = 'AGGGGCT'
%      Flipped = 0;
%    Seq2TDT = flipNregion(Seq2,LeftNTs,RightNTs)
%      Seq2TDT = 'TCGGGGA'
%      Flipped = 1;

function [Nseq, Flipped] = flipNregion(Nseq,varargin)
Flipped = 0;

%Convert Purines AG into X, and Pyrimidine CT into X.
if length(varargin) == 2
    LeftNTs = varargin{1};
    RightNTs = varargin{2};
    
    %Need to insert the "|"
    LeftNT(1:length(LeftNTs)*2-1) = '|';
    LeftNT(1:2:end) = LeftNTs;
    RightNT(1:length(RightNTs)*2-1) = '|';
    RightNT(1:2:end) = RightNTs;   
else
    LeftNT = 'a|g';
    RightNT = 'c|t';
end

LeftCount = length(regexpi(Nseq,LeftNT));
RightCount = length(regexpi(Nseq,RightNT));

if RightCount > LeftCount
    Nseq = seqcomplement(Nseq);
    Flipped = 1;
elseif RightCount == LeftCount
    Flipped = -1;    
end