%findColdHotSpots will find hot spots of AID, and then determine the AA
%mutation frequencies at these hot spots.

%1) Find WRC and GYW hotspots
%2) Calculate the mutation R(i-->b)I , where R is the relative mutations of
%   nt i turning to b (c -- > a, g, t), and I is 1 if this changes the AA,
%   or 0 if this does not change the AA. 
%3) Take sum of R*I to get weight  probability of mutating this residue due
%   to AID SHM.
%4) Compute total R*I within a clonotype, which reflect the number of
%   "opportunities" for mutations.
%5) Compute the actual Mutation of AA due to C mutations
%6) Compare the distbution of C mutations per position 
%   EX:    X axis = Position number where mutatable C was found
%          Y1r = number of exp mutations of the C that can lead to new AA
%          Y1s = number of exp mutations of the C that lead to silent
%                mutations
%          Y2r = number of real mutations of the C that lead to new AA
%          Y2s = number of real mutations that lead to silent mutations

function findColdHotspot(varargin)

[VDJdata, VDJheader] = getPlotVDJdata(varargin);
codonSRMap = getCodonSRMap; %Gets logical arrays indicating which nt in the codong will change the AA
MotifData = getMotifData(VDJdata, VDJheader); %Gets the repertoire motif and rel mutation data
NormMotifData = normalizeMotifData(MotifData);

H = getHeavyHeaderVar(VDJheader);
GrpNum = cell2mat(VDJdata(:, H.GrpNumLoc));
UnqGrpNum = unique(GrpNum)

for y = 1:length(UnqGrpNum)
    GrpIdx = find(GrpNum == UnqGrpNum(y));
    CDR3s = VDJdata{GrpIdx(1), H.CDR3Loc(3)};
    Frame = mod(CDR3s-1, 3) + 1;

    NonsynDist = zeros(1, length(VDJdata{GrpIdx(1), H.RefSeqLoc}));
    SynDist = zeros(1, length(VDJdata{GrpIdx(1), H.RefSeqLoc}));
    for j = 1:length(GrpIdx)
        RefSeq = VDJdata{GrpIdx(j), H.RefSeqLoc};
        Seq = VDJdata{GrpIdx(j), H.SeqLoc};
        
        [~, ALoc, CLoc, GLoc, TLoc] = labelHotSpot(RefSeq);  %Marks the hot spots for mutations
        [SynLoc, NonsynLoc] = labelMutType(RefSeq, Seq, 'Frame', Frame)  %Marsk if these are silent or replacement mutations
        %Binary map of if a -> b leads to silent mutation (0) or replacement (1)
        Igl = zeros(4, length(RefSeq));
        for q = Frame:3:length(RefSeq)-2
            Igl(:, q:q+2) = codonSRMap(Seq(q:q+2));
        end
        
        %Compute expteced AID mutations
        CSynIdx = find(SynLoc & CLoc);
        for k = 1:length(CSynIdx)
            Motif = [RefSeq(CSynIdx(k)-2:CSynIdx(k)) '3'];
            SynExpMut = sum(NormMotifData.Syn.(Motif)' .* ~Igl(:, CSynIdx(k)));
            NonsynExpMut = 1 - SynExpMut;
            SynDist(CSynIdx(k)) = SynDist(CSynIdx(k)) + SynExpMut;
            NonsynDist(CSynIdx(k)) = NonsynDist(CSynIdx(k)) + NonsynExpMut;
        end
        
        GSynIdx = find(SynLoc & GLoc);
        for k = 1:length(GSynIdx)
            Motif = [RefSeq(GSynIdx(k):GSynIdx(k)+2) '1'];
            SynExpMut = sum(NormMotifData.Syn.(Motif)' .* ~Igl(:, GSynIdx(k)));
            NonsynExpMut = 1 - SynExpMut;
            SynDist(GSynIdx(k)) = SynDist(GSynIdx(k)) + SynExpMut;
            NonsynDist(GSynIdx(k)) = NonsynDist(GSynIdx(k)) + NonsynExpMut;
        end
    end
end
