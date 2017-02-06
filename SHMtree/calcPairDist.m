%calcPairDist will compute the pairwise distance matrix and return an
%square matrix, where diagonals are max value or Inf to prevent the
%self-to-self sequence distance from being the lowest distance.
%
%  PairDist = calcPairDist(Seq,DistMode)
%
%  INPUT 
%    Seq: Mx1 cell of trimmed and aligned sequences
%    DistMode: Method for calculating seq-to=seq distance 
%      'ham' for hamming distance
%      'shmham' for hamming distance adjusted based on SHM tendencies. 
%
%  OUTPUT
%    PairDist: MxM matrix of pairwise sequence distance in UINT format.
%    Each row is the parent seq, each col is the child sequence. See note
%    about 'shmham' dist output below.
% 
%  NOTE
%    'X' = wildcard match.
%
%    All sequences should be uppercase.
%
%    When using 'shmham', the PairDist matrix will be double the value to
%    ensure all distance are integers, allowing for low-memory UINT format
%    matricies. The real shmham dist computes a 0.5 correction per
%    mismatched nt that follows SHM tendencies, which would require DOUBLE
%    format matrices that are more memory intensive. This will not affect
%    the lineage tree reconstruction.
%
%  EXAMPLE
%    Seq = {'AAGGTG'; 'ACGCTG'; 'ACTCTG'}
%    [PairDist, HamDist] = calcPairDist(Seq,'shmham')
%    PairDist =
%       255    5   20
%         6  255    3
%        21    3  255
%    HamDist =
%         0    2    3
%         2    0    1
%         3    1    0
function [PairDist, varargout] = calcPairDist(Seq,DistMode,varargin)
%Convert Seq cells to char array
if ischar(Seq)
    SeqT = Seq;
    Seq = cell(size(SeqT,1),1);
    for j = 1:length(Seq)
        Seq{j} = SeqT(j,:);
    end
    clear SeqT;
end

%Determine maximum distance value possible given SeqLength
MaxDist = length(Seq{1})*3; %full mismatch that all agree with SHM tendencies
if MaxDist <= 2^8-1
    Class = 'uint8';
elseif MaxDist <= 2^16-1
    Class = 'uint16';
elseif MaxDist <= 2^32-1
    Class = 'uint32';
end

%Make pairwise difference matrix
PairDist = zeros(length(Seq),Class);
switch lower(DistMode)
    case 'ham' %Hamming distance
        for r = 1:length(Seq)
            Seq1 = Seq{r};
            Xloc1 = Seq1 == 'X';
            for c = 1:r-1
                Seq2 = Seq{c};
                Xloc2 = Seq2 == 'X';
                SeqMatch = (Seq1 == Seq2) | Xloc1 | Xloc2;
                PairDist(r,c) = sum(~SeqMatch);
                PairDist(c,r) = PairDist(r,c);
            end
        end
        HamDist = PairDist; %They are the same
    case 'shmham' %Hamming distance adjusted for SHM tendencies
        HamDist = PairDist;
        for r = 1:length(Seq)
            for c = 1:r-1
                if ~isempty(varargin)
                    [Par2Child, Child2Par, HD] = calcSHMHAMdist(Seq{r},Seq{c},varargin{1});
                else
                    [Par2Child, Child2Par, HD] = calcSHMHAMdist(Seq{r},Seq{c});
                end
                PairDist(r,c) = Par2Child*2; %Multiplying by 2, because calcSHMHAMdist gives you a fraction, and that won't work with uint.
                PairDist(c,r) = Child2Par*2;
                HamDist(r,c) = HD;
                HamDist(c,r) = HD;
            end
        end
end
PairDist(eye(size(PairDist))>0) = Inf; %Set diag to inf to prevent self-pair

if nargout == 2
    varargout{1} = HamDist;
end
