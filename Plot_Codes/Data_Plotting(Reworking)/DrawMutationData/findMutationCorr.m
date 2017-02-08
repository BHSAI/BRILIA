%findMutationCorr will plot the AT mutation% vs GC mutation% to see if
%there is a linear correlation. This would indicate ADAR and AID act same
%time. If you see a vertical / horizontal instead, that would indicate they
%do not act same time. Y axis is ADAR mutation, whereas X axis is AID
%mutation.
%
%We do NOT consider N regions, as there are highly biased towards C/G, and
%thus normalization would be ruined. We thus focus ONLY onlg annotated
%V,D,J regions.

%Output is [AmutCt CmutCt GmutCt Tmutct TotA TotC TotG TotT
%  XmutCt = number of X that were mutated in the RefSeq
%  TotX = number of X that exists in the RefSeq

function Output = findMutationCorr(varargin)
if isempty(varargin)
    [VDJdata, VDJheader, ~, ~] = openSeqData;
else
    VDJdata = varargin{1};
    VDJheader = varargin{2};
end
H = getHeaderVar(VDJheader);

Output = zeros(size(VDJdata,1),8);
for j = 1:size(VDJdata,1)
    %GroupLen = sum((VDJdata{j,H.GrpNumLoc} == GrpNum));
    SamSeq = VDJdata{j,H.SeqLoc}; 
    RefSeq = VDJdata{j,H.RefSeqLoc};    
    if length(RefSeq) ~= length(SamSeq) %Something is wrong, skip
        continue
    end
    
    %Extract only the VDJ nts
    VMDNJ = cell2mat(VDJdata(j,H.LengthLoc));
    VDJloc = [1:VMDNJ(1) sum(VMDNJ(1:2))+1:sum(VMDNJ(1:3)) sum(VMDNJ(1:4))+1:sum(VMDNJ)];
    SamVDJ = SamSeq(VDJloc);
    RefVDJ = RefSeq(VDJloc);
    
    %Extract only the mismatched ones
    MissLoc = RefVDJ ~= SamVDJ;
    RefMiss = upper(RefVDJ(MissLoc));
    
    %Add up total RefSeq A,C,G,T nt count within the VDJ region (no N)
    Act = length(regexpi(RefVDJ,'A'));
    Cct = length(regexpi(RefVDJ,'C'));
    Gct = length(regexpi(RefVDJ,'G'));
    Tct = length(regexpi(RefVDJ,'T'));
    
    %Add up the mutation pairs of highest frequency
    AmutCt = sum(RefMiss == 'A');
    CmutCt = sum(RefMiss == 'C');
    GmutCt = sum(RefMiss == 'G');
    TmutCt = sum(RefMiss == 'T');
    
    Output(j,:) = [AmutCt CmutCt TmutCt GmutCt Act Cct Gct Tct];    
end    
