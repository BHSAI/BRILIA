%flipNregion will analyze a sequence in the N region, and determine if it
%was made from TDT 3' end of the coding or noncoding DNA strand. Returns
%the complement strand IF there are more nucleotides that are assocaited
%with the complement synthesis. The assumption is that elongation of the 3'
%of the coding strand will add G's and A's, whereas elongation of the
%noncoding strand will add C's and T's.
%
%  [Nseq, Flipped] = flipNregion(Nseq)
%
%  [Nseq, Flipped] = flipNregion(Nseq,CodeSideNT,RightNTs)
%
%  INPUT
%    Nseq: nucleotide sequence
%    CodeSideNT: expected nts if TDT acts on + sense strand 3' end (AG)
%    NonCodeSideNT: expected nts if TDT acts on - sense strand 3' end (CT)
%
%  OUTPUT
%    Nseq: flipped (or unflipped) nucleotide sequence
%    Flipped: 1 for flipped, 0 for unflipped status
%
%  EXAMPLE
%    CASE1) Cannot decide how to flip sequence
%      CodeSideNT = 'AG'
%      NonCodeSideNT = 'CT'
%      Seq1 = 'AGGGGCT'
%      [Seq1TDT, Flipped] = flipNregion(Seq1,CodeSideNT,NonCodeSideNT)
%      Seq1TDT = 
%             'AGGGGCT'
%      Flipped = 
%             0
%
%    CASE2) Can decide how to flip sequence
%      CodeSideNT = 'AG'
%      NonCodeSideNT = 'CT'
%      Seq2 = 'AGCCCCT'
%      [Seq2TDT, Flipped] = flipNregion(Seq2,CodeSideNT,NonCodeSideNT)
%      Seq2TDT = 
%             'TCGGGGA'
%      Flipped = 
%             1
%
%  See also calcTDTscore

function [Nseq, Flipped] = flipNregion(Nseq,varargin)
Flipped = 0;

%Convert Purines AG into X, and Pyrimidine CT into X.
if length(varargin) == 2
    CodeSideNT = varargin{1};
    NonCodeSideNT = varargin{2};
    
    %Need to insert the "|"
    CodeNT(1:length(CodeSideNT)*2-1) = '|';
    CodeNT(1:2:end) = CodeSideNT;
    NonCodeNT(1:length(NonCodeSideNT)*2-1) = '|';
    NonCodeNT(1:2:end) = NonCodeSideNT;   
else
    CodeNT = 'a|g';
    NonCodeNT = 'c|t';
end

CodeSideCount = length(regexpi(Nseq,CodeNT));
NonCodeSideCount = length(regexpi(Nseq,NonCodeNT));

if NonCodeSideCount > CodeSideCount
    Nseq = seqcomplement(Nseq);
    Flipped = 1;
elseif NonCodeSideCount == CodeSideCount
    Flipped = -1;    
end
