%findHotSpots will classify the mutation and the surrounding BPs to
%identify if there are any hotspots associated with a point mutation.
%
%  WRCY
%  RGYW 
%  W = A or T
%  R = A or G
%  Y = T or C

%This will look for a NT, then look +- 2 nts afrom the nucleotide and build
%a compositional matrix. It will also produce an 2nd matrix to find
%compositional variation for any C, which is used to determine the degree
%of extra selection power for a hot spot.


function [Amat, Cmat, Gmat, Tmat] = findHotSpots(varargin)
if isempty(varargin)
    [VDJdata, VDJheader, ~, ~] = openSeqData;
else
    VDJdata = varargin{1};
    VDJheader = varargin{2};
end
H = getHeaderVar(VDJheader);

ScanDir = 2;
ScanLen = ScanDir*2+1;

Ahs = zeros(4,ScanLen); %hot spots for mutated ones
Chs = zeros(4,ScanLen); %hot spots for mutated ones
Ghs = zeros(4,ScanLen); %hot spots for mutated ones
Ths = zeros(4,ScanLen); %hot spots for mutated ones
Aas = zeros(4,ScanLen); %all spots available
Cas = zeros(4,ScanLen); %all spots available
Gas = zeros(4,ScanLen); %all spots available
Tas = zeros(4,ScanLen); %all spots available
Ycoor = [1:ScanLen]';

%Determine locations of 1st seq of each group. Want to skip the "ancestral"
%hot spot prediction, because you don't know where they came from.
GrpNum = cell2mat(VDJdata(:,H.GrpNumLoc));
UnqGrpNum = unique(GrpNum);
[~,UnqIdx,~] = intersect(GrpNum,UnqGrpNum); %Tells you where the 1st seq of group starts.

for j = 1:size(VDJdata,1)
    if ~isempty(intersect(j,UnqIdx)) %Skip 1st seq of group, since you don't know how many SHM has occured to germline
        continue
    end
    
    SamSeq = VDJdata{j,H.SeqLoc}; 
    RefSeq = VDJdata{j,H.RefSeqLoc};
    VDJloc = find(SamSeq ~= RefSeq);
    
    %Find only mismatched ones
    VDJloc(VDJloc <= ScanDir ) = []; %Prevent matching left edge nts
    VDJloc(VDJloc >= length(RefSeq)-ScanDir+1) = []; %Prevent matching right edge nts
    for k = 1:length(VDJloc)
        NTref = upper(RefSeq(VDJloc(k)));
        NTsur = nt2int(RefSeq(VDJloc(k)-ScanDir:VDJloc(k)+ScanDir))';
        if max(NTsur) > 4; continue; end %There is an "N" nucleotide.
        Idx = sub2ind([4 ScanLen],NTsur,Ycoor);
        
        switch NTref
            case 'A'
                Ahs(Idx) = Ahs(Idx) + 1;
            case 'C'
                Chs(Idx) = Chs(Idx) + 1;
            case 'G'
                Ghs(Idx) = Ghs(Idx) + 1;
            case 'T'
                Ths(Idx) = Ths(Idx) + 1;
        end
    end
    
    %Find ANY nucleotide
    VDJloc = 1:length(RefSeq);
    VDJloc(VDJloc <= ScanDir ) = [];
    VDJloc(VDJloc >= length(RefSeq)-ScanDir+1) = [];
    for k = 1:length(VDJloc)
        NTref = upper(RefSeq(VDJloc(k)));
        NTsur = nt2int(RefSeq(VDJloc(k)-ScanDir:VDJloc(k)+ScanDir))';
        if max(NTsur) > 4; continue; end %There is an "N" nucleotide.
        Idx = sub2ind([4 ScanLen],NTsur,Ycoor);
        
        switch NTref
            case 'A'
                Aas(Idx) = Aas(Idx) + 1;
            case 'C'
                Cas(Idx) = Cas(Idx) + 1;
            case 'G'
                Gas(Idx) = Gas(Idx) + 1;
            case 'T'
                Tas(Idx) = Tas(Idx) + 1;
        end
    end 
end

Amat = {Ahs Aas};
Cmat = {Chs Cas};
Gmat = {Ghs Gas};
Tmat = {Ths Tas};
