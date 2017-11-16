%calcSHMHAMdist will calculate the adjusted hamming distance that reduces
%the distance if mutations agree with SHM tendencies, and increases the
%distance for consecutively mismatched nts.
%
%Distance = sum(ContigMiss^2) - 0.5sum(FavoredMut) + 0.5sum(UnfavoredMut)
%
%ContigMiss = length of a contiguously mismatched segments
%FavoredMut = number of mutations that agree with SHM tendencies
%UnfavoredMut = number of mutations that disagree with SHM tendencies
%Neutral mutations do not have favored/unfavored distances
%
%The sum(ContigMiss)^2 term was added to increase the distance between
%sequences with many contiguously mismatched nts that are more likely to
%arise from incorrect annotations than SHMs. The idea is by increasing
%distance, they are likely to be cut off from the main tree and reassigned
%a new germline annotation.
%
%  [PCdist, CPdist, HamDist] = calcSHMHAMdist(Seq1, Seq2)
%
%  [PCdist, CPdist, HamDist] = calcSHMHAMdist(Seq1, Seq2, SHMtendency)
%
%  INPUT
%    Seq1: 1st sequence, assumed parent seq
%    Seq2: 2nd sequence, assumed child seq
%    SHMtendency: 4x4 matrix where 1 = favorable mutation, -1 =
%      unfavorable mutations, 0 = no preferrence. Column are ACGT of
%      starting sequence X0, rows are ACGT of ending sequence X1. Default:
%      SHMtendency =    A0  C0  G0  T0
%                    A1  0  -1   1  -1
%                    C1  0   0  -1   0
%                    G1  1  -1   0  -1
%                    T1  1   1  -1   0
%
%  OUTPUT
%    PCdist: distance if Seq1 is parent and Seq2 is child.
%    CPdist: distance if Seq1 is child and Seq2 is parent.
%    HamDist: number of mismatched letters between Seq1 and Seq2.
%
%  NOTE
%    'X' is a wildcard match.
%
%    Low-quality sequence reads can lead to arbitrarily high distances that
%    prevent proper linking of sequences within the same lineage. Consider
%    removing low-quality edges or replacing with wildcard 'X' before
%    processing.
%
%  EXAMPLE
%    Case when C->T and A->G high preference mutations
%      Seq1 = 'ACGTAT' 
%      Seq2 = 'ATGTGT'
%      [PCdist CPdist HamDist] = calcSHMHAMdist(Seq1, Seq2)
%         PCdist =
%             1.0
%         CPdist =
%             1.5
%         HamDist =
%             2
%
%   Case when there is a triple consecutive mismatch  
%     Seq1 = 'ACGCTT' 
%     Seq2 = 'ATTGTT'
%     [PCdist CPdist HamDist] = calcSHMHAMdist(Seq1, Seq2)
%        PCdist %3^2 - 0.5 + 0.5 + 0.5
%            9.5
%        CPdist %3^2 - 0.0 + 0.5 + 0.5
%            10
%        HamDist = 
%            3

function [PCdist, CPdist, HamDist] = calcSHMHAMdist(Seq1, Seq2, varargin)
warning('%s: This is obsolete. Use calcPairDistMEX', mfilename);


%Need Seq1 to be char
if iscell(Seq1)
    Seq1 = Seq1{1};
end

%Want Seq2 to be cell 
if ischar(Seq2) 
    Seq2 = {Seq2};
end

%Use default or custom SHM tendency matrix
if isempty(varargin) || isempty(varargin{1})
    SHMtendency = [...
        0 -1  1 -1;
        0  0 -1  0;
        1 -1  0 -1;
        1  1 -1  0]; %From BRILIA orig paper
else
    SHMtendency = varargin{1};
    if sum(size(SHMtendency) == [4 4]) ~= 2
        error('%s: the SHMtendency matrix should be a 4x4 matrix.', mfilename);
    end
end

PCdist = zeros(length(Seq2), 1);
CPdist = zeros(length(Seq2), 1);
HamDist = zeros(length(Seq2), 1);

for j = 1:length(Seq2)
    MinLen = min(length(Seq1), length(Seq2{j}));
    MissLoc = ~(Seq1(1:MinLen) == Seq2{j}(1:MinLen) | Seq1 == 'X' | Seq2{j} == 'X'); %Oddly, redoing Seq1=='X' is faster than doing it once... 

    if max(MissLoc) == 0
        continue;
    else
        %Find max consec mismatch, and add that.
        Penalty = 0;
        Misses = 0;
        for k = 1:length(MissLoc)-1
            if MissLoc(k) && MissLoc(k+1)
                if Misses == 0
                    Misses = 2;
                else
                    Misses = Misses + 1;
                end
            else
                if Misses > 0
                    Penalty = Penalty + Misses.^2 - Misses; %You are subtracting Misses, since at the end, you don't want to double count the Misses AND Hamming dist.
                    Misses = 0;
                end
            end
        end
        if Misses > 0 %Final Iteration, in case mismatch occurs at end
            Penalty = Misses.^2 + Penalty - Misses;
        end

        RefNT = nt2int(Seq1(MissLoc));
        SamNT = nt2int(Seq2{j}(MissLoc));
        KeepNT = RefNT < 5 & SamNT < 5; %Prevents wildcard matching
        IdxF = sub2ind([4 4], SamNT(KeepNT), RefNT(KeepNT));
        IdxR = sub2ind([4 4], RefNT(KeepNT), SamNT(KeepNT));

        %Calculate the score
        PF = sum(SHMtendency(IdxF));
        PR = sum(SHMtendency(IdxR));
        PCdist(j) = length(IdxF) - 0.5*PF + Penalty;
        CPdist(j) = length(IdxR) - 0.5*PR + Penalty;     
        HamDist(j) = length(IdxF); %Same as hamming distance
    end
end    