%calcPairDist will compute the pairwise distance matrix and return an
%square matrix, where diagonals are max value or Inf to prevent the
%self-to-self sequence distance from being the lowest distance.
%
%  PairDist = calcPairDist(Seq, DistMode)
%
%  PairDist = calcPairDist(Seq, DistMode, Class)
%
%  PairDist = calcPairDist(Seq, DistMode, Class)
%
%  INPUT 
%    Seq: Mx1 cell of trimmed and aligned sequences
%    DistMode: Method for calculating seq-to-seq distance 
%      'ham' for hamming distance
%      'shmham' for hamming distance adjusted based on SHM tendencies. 
%      'identity' for finding distance as (MissedNT) / min([length(SeqA), 
%         length(SeqB)], or sequence identity mismatched of the smaller
%         sequence.
%    Class: Set the class of the resulting PairDist, if you know what it
%      is.
%
%  OUTPUT
%    PairDist: MxM dissimilarity or distance matrix of pairwise sequence
%      comparisons in UINT8, 16, or 32 format. Each row is the parent seq, 
%      each col is the child sequence. See note about 'shmham' dist output
%      below.
% 
%  NOTE
%    'X' = wildcard match.
%
%    All sequences should be uppercase.
%
%    When using 'shmham', the PairDist matrix will be doubled the value to
%    ensure all distance are integers, allowing for low-memory UINT format
%    matrices. The real shmham dist computes a 0.5 correction per
%    mismatched nt that follows SHM tendencies, which would require DOUBLE
%    format matrices that are more memory intensive. This will not affect
%    the lineage tree reconstruction.
%
%  EXAMPLE
%    Seq = {'AAGGTG'; 'ACGCTG'; 'ACTCTG'}
%    [PairDist, HamDist] = calcPairDist(Seq, 'shmham')
%    PairDist =
%       255    5   20
%         6  255    3
%        21    3  255
%    HamDist =
%         0    2    3
%         2    0    1
%         3    1    0
function [PairDist, varargout] = calcPairDist(Seq, DistMode, varargin)
%Convert Seq cells to char array
if ischar(Seq)
    SeqT = Seq;
    Seq = cell(size(SeqT, 1), 1);
    for j = 1:length(Seq)
        Seq{j} = SeqT(j, :);
    end
    clear SeqT;
end
m = length(Seq);

%Ensure DistMode is valid
if ~ismember(DistMode, {'ham', 'shmham', 'identity'})
    error('%s: Unknown DistMode as input.', mfilename);
end

%See if there is a processor number spec in varargin
NumProc = 1;
for j = 1:length(varargin)
    if isnumeric(varargin{j})
        NumProc = varargin{j};
        varargin(j) = [];
    elseif strcmpi(varargin{j}, 'max')
        NumProc = varargin{j};
        varargin(j) = [];
    end
end

%Prevent parfor from triggering auto parallel processing
try
    ps = parallel.Settings;
    ps.Pool.AutoCreate = false;
catch
end

%Setup parallel processing
if NumProc > 1 && ~strcmpi(DistMode, 'ham')
    setParallelProc(NumProc)
end
if NumProc == 1 && ~strcmpi(DistMode, 'ham') && length(Seq) > 100
    disp('Speed up calcPairDist by specifying number of processors > 1')
end

%Determine maximum distance value possible given SeqLength
if isempty(varargin)
    MaxDist = length(Seq{1})^2+length(Seq{1}); %full mismatch that all disagree with SHM tendencies
    if MaxDist <= 2^8-1
        Class = 'uint8';
    elseif MaxDist <= 2^16-1
        Class = 'uint16';
    elseif MaxDist <= 2^32-1
        Class = 'uint32';
    else
        Class = 'double';
    end
else
    Class = lower(varargin{1});
    if ~ismember(Class, {'uint8', 'uint16', 'uint32', 'double'})
        error('%s: Unknown class as input', mfilename);
    end
end

%Make pairwise difference matrix
switch lower(DistMode)
    case 'ham' %Hamming distance (symmetric matrix)
        PairDist = zeros(1, length(Seq)*(length(Seq)-1)/2, Class);
        for c = 1:length(Seq)-1
            S2 = c*(m-1-(c-1)/2);
            S1 = S2 - (m-c) + 1;
            PairDist(S1:S2) = calcHAMdist(Seq(c), Seq(c+1:end));
        end
        PairDist = squareform(PairDist);
        HamDist = PairDist; %They are the same
    
    case 'shmham' %SHM distance (asymmetric matrix)
        PairDist1 = zeros(length(Seq), Class);
        PairDist2 = zeros(length(Seq), Class);
        HamDist1 = zeros(length(Seq), Class);
        HamDist2 = zeros(length(Seq), Class);
        for r = 1:length(Seq)
            %fprintf('%s: %d of %d\n', mfilename, r, length(Seq));
            SeqRow = Seq{r};
            parfor c = 1:r-1
                SeqCol = Seq{c};
                [Par2Child, Child2Par, HD] = calcSHMHAMdist(SeqRow, SeqCol);
                PairDist1(r, c) = Par2Child*2; %Multiplying by 2, because calcSHMHAMdist gives you a fraction, and that won't work with uint.
                PairDist2(c, r) = Child2Par*2;
                HamDist1(r, c) = HD;
                HamDist2(c, r) = HD;
            end
        end
        PairDist = PairDist1 + PairDist2;
        HamDist = HamDist1 + HamDist2;

    case 'identity' %  #missed/lengthofshorterseq
        PairDist = zeros(1, length(Seq)*(length(Seq)-1)/2, 'double');
        HamDist = zeros(1, length(Seq)*(length(Seq)-1)/2, 'uint16');
        for c = 1:length(Seq)-1
            %fprintf('%s: %d of %d\n', mfilename, c, length(Seq));
            Seq1 = Seq{c};
            Seq2 = Seq(c+1:end);
            Dist = zeros(length(Seq2), 1, 'double');
            HDist = zeros(length(Seq2), 1, 'uint16');
            parfor j = 1:length(Seq2)
                Score = alignSeq(Seq1, Seq2{j}, 'Alphabet', 'any');
                MinLen = min([length(Seq1) length(Seq2{j})]);
                Dist(j) = 1 - Score(1)/MinLen;
                HDist(j) = MinLen - Score(1);
            end
            S2 = c*(m-1-(c-1)/2);
            S1 = S2 - (m-c) + 1;
            PairDist(S1:S2) = Dist;
            HamDist(S1:S2) = HDist;
        end
        PairDist = squareform(PairDist);
        HamDist = squareform(HamDist);
end

if nargout == 2
    varargout{1} = HamDist;
end
