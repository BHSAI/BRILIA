%findMutationFreq will look for the nt mutations found between the sequence
%and the germline or predicted sequence. The output is a 4x4 matrix with
%original nt along each column, and mutated nt along each row.
%
%  Hints: Use buildTree before this, if you want parent-child relation.
%  Otherwise, you only calculate difference from germline-child sequence.
%  VDJlen tells the length of the V, D, J segments, cumulative, just in
%  case you want to normalize.
%
%  [Vmat, Dmat, Jmat, VDJlen] = findMutationFreq(VDJdata,NewHeader,Option)
%
%  Option = 'single' will compare mut per sequence
%  Option = 'group' will compare mut per sequence with respect to 1st seq's ref 
function varargout = findMutationFreq(VDJdata,NewHeader,Option)
getHeaderVar;

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
    GrpNum = cell2mat(VDJdata(:,GrpNumLoc));
end
UnqGrpNum = unique(GrpNum);

for y = 1:length(UnqGrpNum)
    IdxLoc = find(UnqGrpNum(y) == GrpNum);
    
    for j = 1:length(IdxLoc)
        %Obtain necessary informations
        Seq = char(VDJdata{IdxLoc(j),SeqLoc});
        RefSeq = char(VDJdata{IdxLoc(1),RefSeqLoc});    
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
        VMDNJ = cell2mat(VDJdata(IdxLoc(j),LengthLoc));
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
