%calcSHMhamdist will calculate an adjusted, hamming distance that is
%reduced for favorable parent-child relation, and increased for otherwise,
%based on the directionality of SHM. 
%
% Score = HamDist + FavorScore + ConsecMiss^2 

%Calculate total probability of SHM by
%multiplying pairwise nucleotide mutation probability P =
%P(N1->M2)xP(N3->M3)x....xP(Nx->Mx), and N is the starting nucleotide, and
%M is the resultant nucleotide. The probability matrix is based off of
%C57BL6 mice data, provided by A collins 2015 paper but processed by
%MAVRIC.
%
% From ancestor ACGT(each col) to child ACGT(each row) mutation matrix
%
% which was determined from non-normalized matrix...
% Praw    = [    0   10395   25057    9981;
%            12809       0    9393   24737;
%            34286   10005       0    7375;
%            15368   30254    5561       0];
% Prawnorm = [
%          0    0.0532    0.1284    0.0511;
%     0.0656         0    0.0481    0.1267;
%     0.1756    0.0512         0    0.0378;
%     0.0787    0.1550    0.0285         0];
% Pmatrix = Prawnorm > 0.25^2 - Prawnorm <= 0.25^2, setting eye(Prawnorm) =
% 0.
%
% EX1: C->T and A->G high preference mutations
%   S1 = 'ACGTAT' 
%   S2 = 'ATGTGT'
%   [PCscore CPscore HamDist] = calcSHMHAMdist(S1,S2)
%      PCscore =
%          1.0
%      CPscore =
%          1.5
%      HamDist =
%          2
% EX2: Triple consec mismatch  
%   S1 = 'ACGCTT' 
%   S2 = 'ATTGTT'
%   [PCscore CPscore HamDist] = calcSHMHAMdist(S1,S2)
%      PCscore = 3^2 - 0.5 + 0.5 + 0.5
%          9.5
%      CPscore = 3^2 - 0.0 + 0.5 + 0.5
%          10
%      HamDist = 
%          3

function [PCscore, CPscore, HamDist] = calcSHMHAMdist(S1,S2,varargin)
%Probability matrix
if isempty(varargin)
    %Default from A collins C57BL6 mice data, processed with MAVRIC
    Pmatrix = [...
        0 -1  1 -1;
        0  0 -1  0;
        1 -1  0 -1;
        1  1 -1  0]; %Turns out, only GA,AG, TC,CT pairs do not help much. But AC,AT does
else
    Pmatrix = varargin{1};
    if sum(size(Pmatrix) == [4 4]) ~= 2
        error('Error: input probability matrix should be 4x4');
    end
end

%Locating the mismathced nts
MissLoc = (S1 ~= S2);

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

if sum(MissLoc) == 0 %Exact match 
    PCscore = 0;
    CPscore = 0;
    HamDist = 0;
else
    RefNT = nt2int(S1(MissLoc));
    SamNT = nt2int(S2(MissLoc));
    KeepNT = RefNT <= 4 & SamNT <= 4; %Prevents wildcard matching
    IdxF = sub2ind([4 4],SamNT(KeepNT),RefNT(KeepNT));
    IdxR = sub2ind([4 4],RefNT(KeepNT),SamNT(KeepNT));

    %Calculate the score
    PF = sum(Pmatrix(IdxF));
    PR = sum(Pmatrix(IdxR));
    PCscore = length(IdxF) - 0.5*PF + PenScore;
    CPscore = length(IdxF) - 0.5*PR + PenScore; 
    HamDist = length(IdxF);
end