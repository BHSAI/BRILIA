%getHotspoutMutData will count the number of hotspot in the CDR3 and FWR3,
%and the observed hotspot mutations in the CDR3 and FWR3, all at the
%nucleotide level.
%
%Hotspot for AID: WRC / GYW
%
%[Hcdr, Hfwr, Ncdr, Nfwr] = getHotspoutMutData(VDJdata, VDJheader)
%
%  INPUT
%
%  OUTPUT
%    Hcdr: number of C/G hotspots in the CDR
%    Hfwr: number of C/G hotspots in the FWR
%    Ncdr: number of mutated hotspot C/G's in the CDR
%    Nfwr: number of mutated hotspot C/G's in the FWR
%
%  See also countHotspotAAMut

function HotspotCount = getHotspoutMutData(varargin)
[VDJdata, VDJheader, FileName, FilePath, varargin] = getPlotVDJdata(varargin{:});
[H, L, Chain] = getAllHeaderVar(VDJheader);

VDJdata = filterVDJdata(VDJdata, VDJheader, 'H-Functional', 'Y');

Hcdr = zeros(size(VDJdata, 1), 1);
Hfwr = zeros(size(VDJdata, 1), 1);
Ncdr = zeros(size(VDJdata, 1), 1);
Nfwr = zeros(size(VDJdata, 1), 1);

GrpNum = cell2mat(VDJdata(:, H.GrpNumLoc));
UnqGrpNum = unique(GrpNum);
GrpData = zeros(length(UnqGrpNum), 5); %GrpSize sum(Hcdr) sum(Hfwr) sum(Ncdr) sum(Nfwr)  per each group
for y = 1:length(UnqGrpNum)
    GrpIdx = find(UnqGrpNum(y) == GrpNum);
    CDR3s = VDJdata{GrpIdx(1), H.CDR3Loc(3)};
    CDR3e = VDJdata{GrpIdx(1), H.CDR3Loc(4)};
    try
        for j = 1:length(GrpIdx)
            RefSeq = VDJdata{j, H.RefSeqLoc};
            Seq = VDJdata{j, H.SeqLoc};
            MutLoc = RefSeq ~= Seq;
            [~, ~, CLoc, GLoc, ~] = labelHotSpot(RefSeq);
            CGLoc = CLoc | GLoc; %All available C or G mutation hotspots
            CGMutLoc = CGLoc & MutLoc; %Hotspot C or G that actually mutated
            Hcdr(GrpIdx(j)) = sum(CGLoc(CDR3s:CDR3e));
            Hfwr(GrpIdx(j)) = sum(CGLoc(1:CDR3s-1));
            Ncdr(GrpIdx(j)) = sum(CGMutLoc(CDR3s:CDR3e));
            Nfwr(GrpIdx(j)) = sum(CGMutLoc(1:CDR3s-1));
        end
        GrpData(y, :) = [length(GrpIdx) sum(Hcdr(GrpIdx)) sum(Hfwr(GrpIdx)) sum(Ncdr(GrpIdx)) sum(Nfwr(GrpIdx))];
    catch
        pause
    end
end

HotspotCount.Hcdr = Hcdr;
HotspotCount.Hfwr = Hfwr;
HotspotCount.Ncdr = Ncdr;
HotspotCount.Nfwr = Nfwr;
HotspotCount.GrpData = GrpData;
