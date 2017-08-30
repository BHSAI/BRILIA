%calcHAMdist will return a matrix with the number of mismatched letters
%between one sequence and one or many other sequence(s). By default, an 'X'
%is a wildcard.
%
%  HamDist = calcHAMdist(Seq1, Seq2)
%
%  HamDist = calcHAMdist(Seq1, Seq2, MatchWild)
%
%  INPUT
%    Seq1: a sequence. UPPERCASE only.
%    Seq2: a sequence or a Mx1 cell of sequences. UPPERCASE only.
%    MatchWild ['y' 'n']: Determines if an 'X' character should be treated
%      as a wildcard or not. Default is MatchWild = 'y'. 
%
%  OUTPUT
%    HamDist: Mx1 matrix of mismatched letters between Seq1 and Seq2.
%
%  NOTE
%    If Seq1 and Seq2 are variable lengths, then HamDist is calculated as
%    number of mismatched letters up to the length of the shorter sequence, 
%    plus the difference in sequence lengths:
%
%  EXAMPLES
%    Seq1 = 'ACXTAT' 
%    Seq2 = 'ATGTGT'
%    HamDist = calcHAMdist(Seq1, Seq2)
%    HamDist = 
%               2
%
%    HamDist = calcHAMdist(Seq1, Seq2, 'n')
%    HamDist = 
%               3
%
%    Seq2 = {'ACXT', 'ATGTGT'};
%    HamDist = calcHAMdist(Seq1, Seq2, 'n')
%    HamDist = 
%               0
%               3
function Dist = calcHAMdist(Seq1, Seq2, varargin)
%Need Seq1 to be char
if iscell(Seq1)
    Seq1 = Seq1{1};
end

%Want Seq2 to be cell 
if ischar(Seq2) 
    Seq2 = {Seq2};
end

%Determine if to allow wild card X matching
if ~isempty(varargin) && ischar(varargin{1}) && ismember(lower(varargin{1}(1)), {'y', 'n'})
    MatchWild = lower(varargin{1}(1));
else
    MatchWild = 'y';
end

%Calc Hamming distances 
Dist = zeros(length(Seq2), 1);
if MatchWild == 'n'
    for j = 1:length(Seq2)
        MinLen = min(length(Seq1), length(Seq2{j}));
        Dist(j) = sum(Seq1(1:MinLen) ~= Seq2{j}(1:MinLen));
    end
else
    XLoc1 = Seq1 == 'X'; %UPPER CASE X ONLY.
    for j = 1:length(Seq2)
        MinLen = min(length(Seq1), length(Seq2{j}));
        XLoc2 = Seq2{j} == 'X';
        Dist(j) = sum(~(Seq1(1:MinLen) == Seq2{j}(1:MinLen) | XLoc1(1:MinLen) | XLoc2(1:MinLen)));
    end
end