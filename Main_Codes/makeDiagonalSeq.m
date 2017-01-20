%makeDiagonalSeq will create a character matrix of Seq2 that spans the
%width of Seq1, such that Diag1 and Diag2 can be compared directly to
%determine alignment. This avoids the iterative sequence checking if one
%were to actually move Seq1 across Seq2 for matches.
%
%  [Diag1, Diag2, GoodIdx] = makeDiagonalSeq(SeqA,SeqB)
%
%  INPUT:
%    SeqA: Sequence of letters
%    SeqB: Sequecne of letters
%
%  OUTPUT
%    Diag1: Repeated char matrix of SeqA
%    Diag2: Char matrix of SeqB but shifted 1 nt per row.
%    GoodIdx: MxN logical matrix indiciating which nt matches should be
%      counted towards the score. Used mainly if you don't want edges to
%      be counted towards the scoring function.
%
%  EXAMPLE
%    SeqA = 'ABCD'
%    SeqB = 'ZXCABCDGH';
%    [Diag1, Diag2, GoodIdx] = makeDiagonalSeq(SeqA,SeqB)
%    Diag1 =
%         ABCD
%         ABCD
%         ABCD
%         ABCD
%         ABCD
%         ABCD
%         ABCD
%         ABCD
%         ABCD
%         ABCD
%         ABCD
%         ABCD
% 
%     Diag2 =
%         ---Z
%         --ZX
%         -ZXC
%         ZXCA
%         XCAB
%         CABC
%         ABCD
%         BCDG
%         CDGH
%         DGH-
%         GH--
%         H---
% 
%     GoodIdx =
%          0     0     0     1
%          0     0     1     1
%          0     1     1     1
%          1     1     1     1
%          1     1     1     1
%          1     1     1     1
%          1     1     1     1
%          1     1     1     1
%          1     1     1     1
%          1     1     1     0
%          1     1     0     0
%          1     0     0     0
%
%  See also convolveSeq, calcAlignScore

function [Diag1, Diag2, GoodIdx] = makeDiagonalSeq(SeqA,SeqB,varargin)
%Determine the number of rows in the diagonal matrix
Tlen = length(SeqA) + length(SeqB) - 1; %Length of the untrimmed diag matrix
if ~isempty(varargin)
    DiagIdx = varargin{1};
else
    DiagIdx = [];
end

%Determine what characters to use for unmatched regions
if ischar(SeqA)
    FillChar = '-';
else
    FillChar = -1;
end

%Initialize diag matrices
Diag1 = repmat(SeqA,Tlen,1);

%Import the diagonalize index here, or create it (slower);
if ~isempty(DiagIdx)
    Diag2Idx = DiagIdx(1:Tlen,end-length(SeqA)+1:end);

%Create the Diag2 matrix from scrap, but it's slower
else 
    Diag2Idx = repmat([-length(SeqA)+2:1],Tlen,1) + repmat([0:Tlen-1]',1,length(SeqA)); %Determine the index for SeqB  %SLOWEST    
end

%Create the Diag2
Diag2 = repmat(FillChar,Tlen,length(SeqA)); %Initalize with fill char
GoodIdx = Diag2Idx > 0 & Diag2Idx <= length(SeqB); %Determine the good index, in binary form.
Diag2(GoodIdx) = SeqB(Diag2Idx(GoodIdx)); %Fill Diag2 with SeqB info.



% Slightly faster, but hard to read
%
%Diag2Idx = ones(Tlen,1)*[-length(SeqA)+2:1] + [0:Tlen-1]'*ones(1,length(SeqA));
% SeqT = [repmat('-',Tlen,1);SeqB'];
% RectSeq = repmat(SeqT,1,Tlen+1);
% RectSeq = RectSeq(:);
% RectSeq(end+1:end+Tlen+1) = '-';
% A = reshape(RectSeq,length(RectSeq)/(Tlen+1),Tlen+1)';
% Diag2 = A(1:Tlen,end-length(SeqB)-length(SeqA)+1:end-length(SeqB));
% GoodIdx = 0;
% Diag2 = Diag2(1:Tlen,1:length(SeqA));
% 
