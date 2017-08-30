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
%         0    5   20
%         6    0    3
%        21    3    0
%    HamDist =
%         0    2    3
%         2    0    1
%         3    1    0
function [PairDist, varargout] = calcPairDist(Seq, DistMode, varargin)
%Convert Seq to a cell array
if ischar(Seq)
    SeqT = Seq;
    Seq = cell(size(SeqT, 1), 1);
    for j = 1:length(Seq)
        Seq{j} = SeqT(j, :);
    end
    clear SeqT;
end
m = length(Seq);

if m <= 1
    PairDist = zeros(m);
    varargout = cell(1, nargout - 1);
    return
end

%Use 16-bit for distance metrics since it's hard to have > 65535 SHMs.
Class = 'uint16';
N = 50; %Number of jobs required to switch from single to parallel core
for j = 1:length(varargin)
    if ischar(varargin{j}) && ismember(lower(varargin{j}), {'uint8', 'uint16', 'uint32', 'double'})
        Class = varargin{j};
    elseif isnumeric(varargin{j})
        N = varargin{j};
    end
end

%Make pairwise difference matrix
switch lower(DistMode)
    case 'ham' %Hamming distance (symmetric matrix)
        PairDist = zeros(1, length(Seq)*(length(Seq)-1)/2, Class);
        if length(Seq) < N %Single core is faster for small jobs
            for r = 1:length(Seq)-1
                S2 = r*(m-1-(r-1)/2);
                S1 = S2 - (m-r) + 1;
                PairDist(S1:S2) = calcHAMdist(Seq{r}, Seq(r+1:end));
            end
        else
            PairDistCell = cell(1:length(Seq) - 1);
            parfor r = 1:length(PairDistCell)
                SeqT = Seq;
                PairDistCell{r} = calcHAMdist(SeqT{r}, SeqT(r+1:end));
            end
            for r = 1:length(Seq)-1
                S2 = r*(m-1-(r-1)/2);
                S1 = S2 - (m-r) + 1;
                PairDist(S1:S2) = PairDistCell{r};
            end
        end
        PairDist = squareform(PairDist);
        HamDist = PairDist; %They are the same
    
    case 'shmham' %SHM distance (asymmetric matrix)
        PairDist = zeros(length(Seq), Class);
        HamDist = zeros(length(Seq), Class);
        if length(Seq) < N %Single core is faster for small jobs
            for r = 1:length(Seq)-1
                [PCD, CPD, HD] = calcSHMHAMdist(Seq{r}, Seq(r+1:end));
                PairDist(r+1:end, r) = CPD*2; %Remember calcSHMHAMdist gives you a half unit, and you want full integers
                PairDist(r, r+1:end) = PCD*2; 
                HamDist(r+1:end, r) = HD;
                HamDist(r, r+1:end) = HD;
            end
        else %Multiple cores are better for larger jobs
            PairDistCell = cell(1, length(Seq) - 1); %Slicable
            parfor r = 1:length(PairDistCell)
                SeqT = Seq; %Just to confirm that Seq should be given to all workers.
                [PCD, CPD, HD] = calcSHMHAMdist(SeqT{r}, SeqT(r+1:end));
                PairDistCell{r} = [PCD*2 CPD*2 HD]; %Remember you want full integer units for PairDist, but true dist has fractions.
            end
            for r = 1:length(Seq) - 1
                PairDist(r, r+1:end) = PairDistCell{r}(:, 1); 
                PairDist(r+1:end, r) = PairDistCell{r}(:, 2); 
                HamDist(r, r+1:end)  = PairDistCell{r}(:, 3);
                HamDist(r+1:end, r)  = PairDistCell{r}(:, 3);
            end
        end
        
    case 'identity' %  #missed/lengthofshorterseq
        PairDist = zeros(1, length(Seq)*(length(Seq)-1)/2, 'double');
        HamDist = zeros(1, length(Seq)*(length(Seq)-1)/2, 'uint16');
        for r = 1:length(Seq)-1
            %fprintf('%s: %d of %d\n', mfilename, c, length(Seq));
            Seq1 = Seq{r};
            Seq2 = Seq(r+1:end);
            Dist = zeros(length(Seq2), 1, 'double');
            HDist = zeros(length(Seq2), 1, 'uint16');
            parfor j = 1:length(Seq2)
                Score = alignSeq(Seq1, Seq2{j}, 'Alphabet', 'any');
                MinLen = min([length(Seq1) length(Seq2{j})]);
                Dist(j) = 1 - Score(1)/MinLen;
                HDist(j) = MinLen - Score(1);
            end
            S2 = r*(m-1-(r-1)/2);
            S1 = S2 - (m-r) + 1;
            PairDist(S1:S2) = Dist;
            HamDist(S1:S2) = HDist;
        end
        PairDist = squareform(PairDist);
        HamDist = squareform(HamDist);
    otherwise
        error('%s: Unrecognized DistMode option.', mfilename);
end

if nargout == 2
    varargout{1} = HamDist;
end
