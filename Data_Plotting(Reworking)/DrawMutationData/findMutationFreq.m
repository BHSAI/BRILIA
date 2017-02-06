%getMutationData will extract the total pairwise mutation data between
%child-parent or child-germline sequences in the VDJdata repertoire data of
%only productive sequences. Use this function to extract the data required
%to use plotting function in this folder.
%
%  MutData = getMutationData
%
%  MutData = getMutationData(VDJdata,VDJheader)
%
%  INPUT
%    VDJdata: MxN cell of all VDJ annotaiton data
%    VDJheader: 1xN cell of all header name for VDJdata
%
%  OUTPUT
%    MutData: structure containing relevant SHM data such. Field are:
%      VmutPar:  V gene 4x4 matrix of X0(child,col) to X1(parent,row)
%      DmutPar:  D gene 4x4 matrix of X0(child,col) to X1(parent,row)
%      JmutPar:  J gene 4x4 matrix of X0(child,col) to X1(parent,row)
%      CDR3mutPar:  CDR3 4x4 matrix of X0(child,col) to X1(parent,row)
%      VmutGerm: V gene 4x4 matrix of X0(child,col) to X1(germline,row)
%      DmutGerm: D gene 4x4 matrix of X0(child,col) to X1(germline,row)
%      JmutGerm: J gene 4x4 matrix of X0(child,col) to X1(germline,row)
%      CDR3mutGerm:  CDR3 4x4 matrix of X0(child,col) to X1(germline,row)
%
%  NOTE
%    This will calculated SHM for ALL sequences in the VDJdata file, except
%    those with unresolved parent sequences. If you want to do this only
%    for productive sequence, filter VDJdata prior to running this
%    function.
%
%    The 4x4 matrix contains the following data for X0 -> X1 mutations,
%    where N is the total count of each mutation occurence. Divide N by the
%    approriate totNT variable to normalize to mutation rate (e.g. mutation
%    chance of A->C per bp). 
%       A0  C0  G0  T0
%    A1  0   N   N   N
%    C1  N   0   N   N
%    G1  N   N   0   N
%    T1  N   N   N   0

function MutData = findMutationFreq(varargin)
%See if user gave VDJdata
if isempty(varargin) || (~isempty(varargin) && isempty(varargin{1})) %Need to find file
    [VDJdata,VDJheader] = openSeqData;
elseif ~isempty(varargin)
    if ischar(varargin{1}) %Filename was given
        [VDJdata,VDJheader] = openSeqData(varargin{1});
    elseif length(varargin) == 2 && iscell(varargin{1}) && iscell(varargin{2}) %VDJdata and VDJheader was given
        VDJdata = varargin{1};
        VDJheader = varargin{2};
    end
else
    error('getMutationData: Check the inputs');
end
H = getHeaderVar(VDJheader);

%Initialize the temporary 3D matrix for storing N number of pairwise muts
ParMutMat = zeros(4,4,5);
GermMutMat = zeros(4,4,5);

%Extract sequence grouping information
GrpNum = cell2mat(VDJdata(:,H.GrpNumLoc));
UnqGrpNum = unique(GrpNum);
for y = 1:length(UnqGrpNum)
    IdxLoc = find(UnqGrpNum(y) == GrpNum);
    GermSeq = VDJdata{IdxLoc(1),H.RefSeqLoc};
    DelIdx3 = regexpi(GermSeq,'[^ACGTU]');
    VMDNJ = cell2mat(VDJdata(IdxLoc(1),H.LengthLoc));
    
    for k = 1:2     %evaluate to ref seq then germ seq
        for j = 1:length(IdxLoc)    %Process each seq in group
            Seq = VDJdata{IdxLoc(j),H.SeqLoc};           
            DelIdx1 = regexpi(Seq,'[^ACGTU]');
            if sum(VMDNJ) ~= length(Seq) %Another bug...
                warning('VMDNJ lengths do not match with seq length in GrpNum %d',UnqGrpNum(y));
                continue 
            elseif isempty(VMDNJ) %Incomplete annotation
                warning('VMDNJ lengths are not resolved in GrpNum %d',UnqGrpNum(y));
                continue
            end

            if k == 1
                RefSeq = VDJdata{IdxLoc(j),H.RefSeqLoc}; %Immediate parent seq
                DelIdx2 = regexpi(RefSeq,'[^ACGTU]');
            else
                RefSeq = GermSeq; %Germ seq
                DelIdx2 = DelIdx3;
            end
            if length(RefSeq) ~= length(Seq)%Just in case there's a bug
                warning('Ref seq length mismatch error in GrpNum %d',UnqGrpNum(y));
                continue
            end

            %Get rid of ambiguous sequences mismatches
            MissIdxb = Seq ~= RefSeq; %binary index of mismatched nt
            MissIdxb([DelIdx1 DelIdx2]) = 0; 
            MissIdx = find(MissIdxb == 1);
            if isempty(MissIdx); continue; end

            %Determine the mutation for the various segments
            VDJmat = zeros(4,4,5); %4x4 matrix for the 5 VMDNJ lengths
            for t = 1:size(VDJmat,3)
                if VMDNJ(t) == 0; continue; end
                if t == 1
                    S1 = 1;
                    S2 = VMDNJ(1);
                else
                    S1 = sum(VMDNJ(1:t-1)) + 1;
                    S2 = sum(VMDNJ(1:t));
                end
                EvalMissIdx = MissIdx(MissIdx>=S1 & MissIdx <= S2);
                
                if ~isempty(EvalMissIdx)
                    X0 = nt2int(Seq(:,EvalMissIdx)); %Tells column locations
                    X1 = nt2int(RefSeq(:,EvalMissIdx)); %Tells row location
                    for q = 1:length(X0)
                        VDJmat(X1(q),X0(q),t) = VDJmat(X1(q),X0(q),t) + 1;
                    end
                end
            end
            
            if k == 1 %Update the structures for child to parent comparison
                ParMutMat = ParMutMat + VDJmat;
            else %Update structures for child to germline comparison
                GermMutMat = GermMutMat + VDJmat;
            end
        end
    end
end

%Initialize the output structure
MutData.VmutPar = ParMutMat(:,1);
MutData.MmutPar = ParMutMat(:,2);
MutData.DmutPar = ParMutMat(:,3);
MutData.NmutPar = ParMutMat(:,4);
MutData.JmutPar = ParMutMat(:,5);
MutData.VmutGerm = GermMutMat(:,1);
MutData.MmutGerm = GermMutMat(:,2);
MutData.DmutGerm = GermMutMat(:,3);
MutData.NmutGerm = GermMutMat(:,4);
MutData.JmutGerm = GermMutMat(:,5);
