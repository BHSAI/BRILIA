%collectMutFreqFromCDR3 will collect the SHM frequencies as determined from
%the start and end of the CDR3 nuceleotide. Distance N for V side means 1
%nt away from the 1st nt of the 105th residues codon. Distance M for J
%side means 1 nt away from the last nt of the 117th residue codon. Distance
%are computed PER each unique V and J gene.
%
function MutDistData = collectMutDistFromCDR3(varargin)
[VDJdata, VDJheader, ~, ~, varargin] = getPlotVDJdata(varargin{:});
P = inputParser;
addParameter(P, 'Chain', 'H', @(x) ismember({upper(x)}, {'H', 'L'}));
parse(P, varargin{:})
Chain = P.Results.Chain;

%Determine which chain datat to extract
[H, L, ChainT] = getAllHeaderVar(VDJheader);
if ~contains(ChainT, Chain)
    warning('%s: No data for this Ig chain - %s', mfilename, Chain);
    return
end
if strcmpi(Chain, 'H')
    B = H;
else
    B = L;
end

%Determine unique V names
Vnames = VDJdata(:, B.GeneNameLoc(1));
UnqVnames = unique(Vnames);
for j = 1:length(UnqVnames)
    Vsplit = regexp(UnqVnames{j}, '\|', 'split');
    if length(Vsplit) > 1
        UnqVnames{j} = Vsplit{1}; %Take first one only
        %UnqVnames = cat(1, UnqVnames, Vsplit(2:end)');
    end
end
UnqVnames = unique(UnqVnames);

%Determine unique J names
Jnames = VDJdata(:, B.GeneNameLoc(end));
UnqJnames = unique(Jnames);
for j = 1:length(UnqJnames)
    Jsplit = regexp(UnqJnames{j}, '\|', 'split');
    if length(Jsplit) > 1
        UnqJnames{j} = Jsplit{1}; %Take first one only
        %UnqJnames = cat(1, UnqJnames, Jsplit(2:end)');
    end
end
UnqJnames = unique(UnqJnames);

%Set arbitrary distance of 400 nts left and right of the CDR region.
Vmut = zeros(length(UnqVnames), 400);
Jmut = zeros(length(UnqJnames), 400);

VidxMap = containers.Map;
for j = 1:length(UnqVnames)
    VidxMap(UnqVnames{j}) = j;
end

JidxMap = containers.Map;
for j = 1:length(UnqJnames)
    JidxMap(UnqJnames{j}) = j;
end

GrpNum = cell2mat(VDJdata(:, B.GrpNumLoc));
UnqGrpNum = unique(GrpNum);
for y = 1:length(UnqGrpNum)
    GrpIdx = find(UnqGrpNum(y) == GrpNum);
    CDR3s = VDJdata{GrpIdx(1), B.CDR3Loc(3)};
    CDR3e = VDJdata{GrpIdx(1), B.CDR3Loc(4)};
    Vname = VDJdata{GrpIdx(1), B.GeneNameLoc(1)};
    Vname = regexp(Vname, '\|', 'split');
    Vname = Vname{1};
    Vidx = VidxMap(Vname);
    
    Jname = VDJdata{GrpIdx(1), B.GeneNameLoc(end)};
    Jname = regexp(Jname, '\|', 'split');
    Jname = Jname{1};
    Jidx = JidxMap(Jname);
    for g = 1:length(GrpIdx)
        RefSeq = VDJdata{GrpIdx(g), B.RefSeqLoc};
        Seq = VDJdata{GrpIdx(g), B.SeqLoc};
        MutLoc = ~(RefSeq == Seq | RefSeq == 'N' | Seq == 'N');
        VmutLoc = fliplr(MutLoc(1:CDR3s-1)); %Want left side to be pos 1, so flip
        JmutLoc = MutLoc(CDR3e+1:end);
        
        Vmut(Vidx, 1:length(VmutLoc)) = Vmut(Vidx, 1:length(VmutLoc)) + VmutLoc;
        Jmut(Jidx, 1:length(JmutLoc)) = Jmut(Jidx, 1:length(JmutLoc)) + JmutLoc;
    end
end

%Trim off exccess Vmut and Jmut
VmutSum = sum(Vmut, 1);
VcutLoc = size(Vmut, 2);
while VmutSum(VcutLoc) == 0
    VcutLoc = VcutLoc - 1;
    if VcutLoc < 1
        break
    end
end
Vmut(:, VcutLoc+1:end) = [];

%Trim off exccess Vmut and Jmut
JmutSum = sum(Jmut, 1);
JcutLoc = size(Jmut, 2);
while JmutSum(JcutLoc) == 0
    JcutLoc = JcutLoc - 1;
    if JcutLoc < 1
        break
    end
end
Jmut(:, JcutLoc+1:end) = [];

%Package the outputs
MutDistData.V = cat(2, UnqVnames, num2cell(Vmut));
MutDistData.J = cat(2, UnqJnames, num2cell(Jmut));
