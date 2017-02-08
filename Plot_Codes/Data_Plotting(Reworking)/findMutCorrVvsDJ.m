%findMutCorrVvsDJ will look at X -> Y nucleotide mutation frequencies
%across V and D+J segments, and then return the % of X that mutated to
%A,C,G,T, etc. Returns a 3x8 matrix, with 
%1st row = [Av Cv Gv Tv Adj Cdj Gdj Tdj] aveerage mutation frequencies for the entire data,
%2nd row = standard deviations.
%The data ONLY includes mutated nts, and will not be diluted by 0
%mutations. Also, if a segment lacks a nucleotide, EX: Dseg = 'GGGG', then
%it will assign a negative frequency for mark for exclusion in calculations of Avg and Std mut rates.

%EX: Vseg = 'TGTAGG', DJseg = 'GGGTGG'
%CorrMat = findMutCorrVvsDJ(VDJdata,VDJheader
%CorrMat = []


function varargout = findMutationFreq(VDJdata,VDJheader,Option)
H = getHeaderVar(VDJheader);

%Matrix for: Ac Ag At Ca Ct Cg Ga Gc Gt Ta Tc Tg Atot Ctot Gtot Ttot
CountMatV = zeros(size(VDJdata,1),16+4); 
CountMatD = zeros(size(VDJdata,1),16+4);
CountMatJ = zeros(size(VDJdata,1),16+4);

for j = 1:size(VDJdata,1)
    %Find the mismatch locations with respect to par-child relation
    VMDNJ = cell2mat(VDJdata(j,H.LengthLoc));
    RefSeq = VDJdata{j,H.RefSeqLoc};
    CurSeq = VDJdata{j,H.SeqLoc};
    MissLoc = find(CurSeq ~= RefSeq);
    
    %Make sure ambiguous nts are excluded from mutation analysis
    AmbLoc0 = regexpi(RefSeq,'[^ACGTU]');
    AmbLoc1 = regexpi(CurSeq,'[^ACGTU]');
    AmbLoc = unique([AmbLoc0 AmbLoc1]);
    MissLoc(AmbLoc) = 0; 
    
    %Start counting pairwise mutations
    H.VmutLoc = MissLoc(MissLoc <= VMDNJ(1));
    H.DmutLoc = MissLoc(MissLoc > sum(VMDNJ(1:2)) & MissLoc <= sum(VMDNJ(1:3)));
    H.JmutLoc = MissLoc(MissLoc > sum(VMDNJ(1:4)));
        
    %Fill in the mutation 4x4 matrix, where col is par, row is child
    Vmut0 = nt2int(RefSeq(H.VmutLoc));
    Vmut1 = nt2int(CurSeq(H.VmutLoc));
    Vmat = zeros(4);
    for x = 1:length(Vmut0)
        Vmat(Vmut1(x),Vmut0(x)) = Vmat(Vmut1(x),Vmut0(x)) + 1;
    end
    
    Dmut0 = nt2int(RefSeq(H.DmutLoc));
    Dmut1 = nt2int(CurSeq(H.DmutLoc));
    Dmat = zeros(4);
    for x = 1:length(Dmut0)
        Dmat(Dmut1(x),Dmut0(x)) = Dmat(Dmut1(x),Dmut0(x)) + 1;
    end
    
    Jmut0 = nt2int(RefSeq(H.JmutLoc));
    Jmut1 = nt2int(CurSeq(H.JmutLoc));
    Jmat = zeros(4);
    for x = 1:length(Jmut0)
        Jmat(Jmut1(x),Jmut0(x)) = Jmat(Jmut1(x),Jmut0(x)) + 1;
    end
    
    %Get the base compositoin of ref sequences
    Vbase = cell2mat(struct2cell(basecount(RefSeq(1:VMDNJ(1)))));
    Dbase = cell2mat(struct2cell(basecount(RefSeq(sum(VMDNJ(1:2))+1:sum(VMDNJ(1:3))))));
    Jbase = cell2mat(struct2cell(basecount(RefSeq(sum(VMDNJ(1:4))+1:sum(VMDNJ(1:5))))));
    
    %Fill in the tables
    CountMatV(j,:) = [Vmat(:); Vbase];
    CountMatD(j,:) = [Dmat(:); Dbase];
    CountMatJ(j,:) = [Jmat(:); Jbase];
end

%Here on, we are doing this for the intented V vs DJ comparisons
CountMatDJ = CountMatD+CountMatJ;

%Look for no X compositions, and mark for deletion
DelV = (CountMatV(:,17) == 0) | (CountMatV(:,18) == 0) | (CountMatV(:,19) == 0) | (CountMatV(:,20) == 0);
DelDJ = (CountMatDJ(:,17) == 0) | (CountMatDJ(:,18) == 0) | (CountMatDJ(:,19) == 0) | (CountMatDJ(:,20) == 0);
DelVDJ = DelV | DelDJ;
CountMatV(DelVDJ,:) = [];
CountMatDJ(DelVDJ,:) = [];

%Now normalize is mutation by number of corresponding X ref nts
CountMatV(:,1:4) = CountMatV(:,1:4)./repmat(CountMatV(:,17),1,4);
CountMatV(:,5:8) = CountMatV(:,5:8)./repmat(CountMatV(:,18),1,4);
CountMatV(:,9:12) = CountMatV(:,9:12)./repmat(CountMatV(:,19),1,4);
CountMatV(:,13:16) = CountMatV(:,13:16)./repmat(CountMatV(:,20),1,4);

CountMatDJ(:,1:4) = CountMatDJ(:,1:4)./repmat(CountMatDJ(:,17),1,4);
CountMatDJ(:,5:8) = CountMatDJ(:,5:8)./repmat(CountMatDJ(:,18),1,4);
CountMatDJ(:,9:12) = CountMatDJ(:,9:12)./repmat(CountMatDJ(:,19),1,4);
CountMatDJ(:,13:16) = CountMatDJ(:,13:16)./repmat(CountMatDJ(:,20),1,4);

%Now delete the aa cc gg tt mutations.
CountMatV(:,[1 6 11 16]) = [];
CountMatDJ(:,[1 6 11 16]) = [];

VtotRate = zeros(2,12);
DJtotRate = zeros(2,12);
for k = 1:12
    KeepThis = (CountMatV(:,k) > 0);
    VtotRate(1,k) = mean(CountMatV(KeepThis,k));
    VtotRate(2,k) = std(CountMatV(KeepThis,k));

    KeepThis = (CountMatDJ(:,k) > 0);
    DJtotRate(1,k) = mean(CountMatDJ(KeepThis,k));
    DJtotRate(2,k) = std(CountMatDJ(KeepThis,k));
end

%Now normalize each to its own respectiv A->X, C->X, G->X, T->X.
VtotRate(1,1:3) = VtotRate(1,1:3)/sum(VtotRate(1,1:3));
VtotRate(1,4:6) = VtotRate(1,4:6)/sum(VtotRate(1,4:6));
VtotRate(1,7:9) = VtotRate(1,7:9)/sum(VtotRate(1,7:9));

DJtotRate(1,1:3) = DJtotRate(1,1:3)/sum(DJtotRate(1,1:3));
DJtotRate(1,4:6) = DJtotRate(1,4:6)/sum(DJtotRate(1,4:6));
DJtotRate(1,7:9) = DJtotRate(1,7:9)/sum(DJtotRate(1,7:9));

scatter(VtotRate(1,:),DJtotRate(1,:))






VDJmat = zeros(4,4,3);
VDJlen = zeros(1,3);

VmutCt = zeros(size(VDJdata,1),1);
MmutCt = zeros(size(VDJdata,1),1);
DmutCt = zeros(size(VDJdata,1),1);
NmutCt = zeros(size(VDJdata,1),1);
JmutCt = zeros(size(VDJdata,1),1);

if strcmpi(Option,'single')
    GrpNum = [1:size(VDJdata,1)]';
else
    GrpNum = cell2mat(VDJdata(:,H.GrpNumLoc));
end
UnqGrpNum = unique(GrpNum);

for y = 1:length(UnqGrpNum)
    IdxLoc = find(UnqGrpNum(y) == GrpNum);
    
    for j = 1:length(IdxLoc)
        %Obtain necessary informations
        Seq = char(VDJdata{IdxLoc(j),H.SeqLoc});
        RefSeq = char(VDJdata{IdxLoc(1),H.RefSeqLoc});    
        if length(RefSeq) ~= length(Seq)
            continue
        end

        %Look for those non-nucleotide char, like X or N
        MismatchLoc = Seq ~= RefSeq;
        DelIdx1 = regexpi(Seq,'[^ACGTU]');
        DelIdx2 = regexpi(RefSeq,'[^ACGTU]');
        MismatchLoc([DelIdx1 DelIdx2]) = 0; %Get rid of ambiguous sequences;
        MissLoc = find(MismatchLoc == 1);

        %Determine the V, D, J locations
        VMDNJ = cell2mat(VDJdata(IdxLoc(j),H.LengthLoc));
        if sum(VMDNJ) ~= length(Seq)
            continue
        elseif isempty(VMDNJ)
            continue
        end
        VDJlen = VDJlen + VMDNJ(1:2:5); %Add up the V,D,J gene lengths

        %Determine the mutation for the various segments
        VDJmissLoc = cell(1,3);
        VDJmissLoc{1} = MissLoc(MissLoc <= VMDNJ(1));
        VDJmissLoc{2} = MissLoc(MissLoc > sum(VMDNJ(1:2)) & MissLoc <= sum(VMDNJ(1:3)));
        VDJmissLoc{3} = MissLoc(MissLoc > sum(VMDNJ(1:4)));
        MmissLoc = MissLoc(MissLoc > VMDNJ(1) & MissLoc <= sum(VMDNJ(1:2)));
        NmissLoc = MissLoc(MissLoc > sum(VMDNJ(1:3)) & MissLoc <= sum(VMDNJ(1:4)));
        for k = 1:3
            if ~isempty(VDJmissLoc{k})
                SeqNT = nt2int(Seq(:,VDJmissLoc{k}));
                RefNT = nt2int(RefSeq(:,VDJmissLoc{k}));
                for q = 1:length(SeqNT)
                    VDJmat(SeqNT(q),RefNT(q),k) = VDJmat(SeqNT(q),RefNT(q),k) + 1;
                end
            end
        end

        %Sum up the total mutations be VMDNJ segment
        VmutCt(IdxLoc(j)) = length(VDJmissLoc{1});
        MmutCt(IdxLoc(j)) = length(MmissLoc);
        DmutCt(IdxLoc(j)) = length(VDJmissLoc{2});
        NmutCt(IdxLoc(j)) = length(NmissLoc);
        JmutCt(IdxLoc(j)) = length(VDJmissLoc{3});    
    end
end
varargout{1} = VDJmat(:,:,1);
varargout{2} = VDJmat(:,:,2);
varargout{3} = VDJmat(:,:,3);
varargout{4} = VDJlen;
varargout{5} = [VmutCt MmutCt DmutCt NmutCt JmutCt];
