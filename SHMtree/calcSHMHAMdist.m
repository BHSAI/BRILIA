%calcSHMHAMdist will calculate the adjusted hamming distance that has
%reduced distance if mutations agree with SHM tendencies, and increased
%distance for consecutively mismatched nts.
%
%Distance = sum(ContigMiss^2) - FavorableMut , 
%where ContigMiss is the number of contiguously mismatched nts, and
%FavorableMut is the adjustment for mutations that agree or disagree with
%SHM tendencies. Each favorable mutation decrease Distance by 0.5,
%unfavorable increase Distance by 0.5, and others do not change Distance.
%
%  [PCdist, CPdist, HamDist] = calcSHMHAMdist(Seq1,Seq2)
%
%  [PCdist, CPdist, HamDist] = calcSHMHAMdist(Seq1,Seq2,SHMtendency)
%
%  INPUT
%    Seq1 = 1st sequence, assumed parent seq
%    Seq2 = 2nd sequence, assumed child seq
%    SHMtendency = 4x4 matrix where 1 = favorable mutation, -1 =
%      unfavorable mutations, 0 = no preferrence. Column are ACGT of
%      starting sequence X0, rows are ACGT of ending sequence X1. Default:
%      SHMtendency =    A0  C0  G0  T0
%                    A1  0  -1   1  -1
%                    C1  0   0  -1   0
%                    G1  1  -1   0  -1
%                    T1  1   1  -1   0
%
%  OUTPUT
%    PCdist = distance if Seq1 is parent and Seq2 is child.
%    CPdist = distance if Seq2 is parent and Seq1 is child.
%    HamDist = number of mismatched letters between Seq1 and Seq2.
%
%  NOTE
%    'X' is a wildcard match.
%
%  EXAMPLE
%    Case when C->T and A->G high preference mutations
%      Seq1 = 'ACGTAT' 
%      Seq2 = 'ATGTGT'
%      [PCdist CPdist HamDist] = calcSHMHAMdist(Seq1,Seq2)
%         PCdist =
%             1.0
%         CPdist =
%             1.5
%         HamDist =
%             2
%   Case when there is a triple consec mismatch  
%     Seq1 = 'ACGCTT' 
%     Seq2 = 'ATTGTT'
%     [PCdist CPdist HamDist] = calcSHMHAMdist(Seq1,Seq2)
%        PCdist %3^2 - 0.5 + 0.5 + 0.5
%            9.5
%        CPdist %3^2 - 0.0 + 0.5 + 0.5
%            10
%        HamDist = 
%            3

function [PCdist, CPdist, HamDist] = calcSHMHAMdist(Seq1,Seq2,varargin)
%Probability matrix
if isempty(varargin)
    %Values obtained from BRILIA paper
    SHMtendency = [...
        0 -1  1 -1;
        0  0 -1  0;
        1 -1  0 -1;
        1  1 -1  0];
else
    SHMtendency = varargin{1};
    if sum(size(SHMtendency) == [4 4]) ~= 2
        error('Error: input probability matrix should be 4x4');
    end
end

%Locating the mismathced nts
MissLoc = ~((Seq1 == Seq2) | (Seq1 == 'X') | (Seq2 == 'X'));

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
    PCdist = 0;
    CPdist = 0;
    HamDist = 0;
else
    RefNT = nt2int(Seq1(MissLoc));
    SamNT = nt2int(Seq2(MissLoc));
    KeepNT = RefNT <= 4 & SamNT <= 4; %Prevents wildcard matching
    IdxF = sub2ind([4 4],SamNT(KeepNT),RefNT(KeepNT));
    IdxR = sub2ind([4 4],RefNT(KeepNT),SamNT(KeepNT));

    %Calculate the score
    PF = sum(SHMtendency(IdxF));
    PR = sum(SHMtendency(IdxR));
    PCdist = length(IdxF) - 0.5*PF + PenScore;
    CPdist = length(IdxF) - 0.5*PR + PenScore; 
    HamDist = length(IdxF);
end