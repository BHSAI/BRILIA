%calcPairDist(Seq,DistMode) will compute the pairwise distance matrix and
%return an NxN matrix, where diagonals are Inf. 
%
%  PairDist = calcPairDist(Seq,DistMode) where Seq is a cell of aligned,
%  tirmmed sequences.
%    DistMode = 'ham' for hamming distance
%    DistMode = 'shm' for SHM distance
%    DistMode = 'shmham' for hamming distance adjusted based on
%    likelihood of Seq1 --> Seq2 or vice versa.

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
SeqLen = size(Seq{1},2)*3;
if SeqLen <= 2^8-1
    Class = 'uint8';
elseif SeqLen <= 2^16-1
    Class = 'uint16';
elseif SeqLen <= 2^32-1
    Class = 'uint32';
end

%Make pairwise difference matrix
PairDist = zeros(size(Seq,1),Class);
switch lower(DistMode)
    case 'ham' %Hamming distance
        for r = 1:size(Seq,1)
            for c = 1:r-1
                PairDist(r,c) = sum(Seq{r} ~= Seq{c});
                PairDist(c,r) = PairDist(r,c);
            end
        end
        HamDist = PairDist;
    case 'hampen' %Hamming distance with penalty for consec mismatches
        for r = 1:size(Seq,1)
            for c = 1:r-1
                MissLoc = Seq{r} ~= Seq{c};
                %Find max consec mismatch, and add that.
                PenScore = 0; %Penalty score
                MissCount = 0;
                for k = 1:length(MissLoc)-1
                    if MissLoc(k) && MissLoc(k+1)
                        if MissCount == 0
                            MissCount = 2;
                        else
                            MissCount = MissCount+1;
                        end
                    else
                        if MissCount > 0
                            PenScore = MissCount.^2 + PenScore - MissCount; %You are subtracting MissCount, since at the end, you don't want to double count the Miss Count AND Ham Dist.
                            MissCount = 0;
                        end
                    end
                end
                if MissCount > 0 %Final Iteration, in case mismatch occurs at end
                    PenScore = MissCount.^2 + PenScore - MissCount;
                end
                
                %Final Score, add the hamming distance
                FinalScore = PenScore + sum(MissLoc);
                
                PairDist(r,c) = FinalScore;
                PairDist(c,r) = FinalScore;
            end
        end
        HamDist = PairDist;
%     case 'shm' %SHM distance
%         HamDist = zeros(size(Seq,1),Class);
%         for r = 1:size(Seq,1)
%             for c = 1:size(Seq,1)
%                 if r ~= c
%                     [PairDist(r,c), ~, HamDist(r,c)] = calcSHMdist(Seq{r},Seq{c});
%                 end
%             end
%         end
    case 'shmham' %SHM distance, averaged over HAM distance
        HamDist = zeros(size(Seq,1),Class);
        for r = 1:size(Seq,1)
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