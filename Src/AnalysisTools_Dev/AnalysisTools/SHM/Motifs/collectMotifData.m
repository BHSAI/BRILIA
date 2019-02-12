%collectMotifData will collect from VDJdata all the mutations that is
%occuring around the various 2- or 3- nucleotid motifs, storing also the
%relative mutation frequencies from X -> Y. This is computed per
%parent-child relationship.
function varargout = collectMotifData(varargin)
P = inputParser;
addParameter(P, 'Chain', '', @(x) ismember({upper(x)}, {'H', 'L', 'HL', ''}));
addParameter(P, 'MutType', 'All', @(x) ismember({upper(x)}, {'All', 'S', 'N'}));
[VDJdata, VDJheader, ~, ~, varargin] = getPlotVDJdata(varargin{:});
[Ps, Pu, ReturnThis, ExpPs, ExpPu] = parseInput(P, varargin{:});
if ReturnThis
   varargout = {Ps, Pu, ExpPs, ExpPu};
   return;
end
Chain = upper(P.Results.Chain);
MutType = P.Results.MutType;

%Create the motif list, given as NNN#
CodonList = loadCodonList;
codonMap = containers.Map(CodonList, 1:length(CodonList));
RelMutCell = repmat({zeros(1, 4)}, length(CodonList), 3); %Stores the X -> Y (1x4) mutations matrix per codon per 3 positions

%Format into a Mx5 cell with TriNucCode, Amut, Cmut, Gmut, Tmut
TriNuc = repmat(CodonList, 1, 3);
for j = 1:length(TriNuc)
    for k = 1:3
        TriNuc{j, k} = [TriNuc{j, k} num2str(k)];
    end
end
MotifData.ColHeader = {'Pos1', 'Pos2', 'Pos3'};
MotifData.RowHeader = CodonList;
MotifData.Data = [];

%Extract the index of necessary data
Map = getVDJmapper(VDJheader);
if isempty(Chain); Chain = Map.Chain; end
[SeqIdx, RefSeqIdx, CDR3SIdx] = deal(zeros(1, 2));
for c = 1:numel(Chain)
    C = lower(Chain(c));
    SeqIdx(c)    = Map.([C 'Seq']);
    RefSeqIdx(c) = Map.([C 'RefSeq']);
    CDR3SIdx(c)  = Map.([C 'CDR3'])(3);
end
SeqIdx    = nonzeros(SeqIdx);
RefSeqIdx = nonzeros(RefSeqIdx);
CDR3SIdx  = nonzeros(CDR3SIdx);
if isempty(SeqIdx) || isempty(RefSeqIdx) || isempty(CDR3SIdx)
    error('%s: Invalid data file. Cannot find correct indices', mfilename);
end

%Create the mutation matrix
for j = 1:size(VDJdata,1)
    MutInfo = [];
    for b = 1:length(RefSeqIdx)
        RefSeq = VDJdata{j, RefSeqIdx(b)};
        Seq = VDJdata{j, SeqIdx(b)};
        CDR3s = VDJdata{j, CDR3SIdx(b)};
        Frame = mod(CDR3s-1, 3) + 1;
        MutInfo = cat(2, MutInfo, getMutInfo(RefSeq, Seq, Frame)); %No frame needed?
    end
    if strcmpi(MutType, 'all')
        for r = 1:length(MutInfo)
            MutMotif = MutInfo(r).Motif;
            IntAfter = nt2int(MutInfo(r).NT(end));
            for k = 1:3
                try
                    CodonLoc = codonMap(MutMotif(5-k-1:5-k+1));
                catch
                    continue;
                end
                RelMutCell{CodonLoc, k}(IntAfter) = RelMutCell{CodonLoc, k}(IntAfter) + 1; 
            end
        end
    else
        for r = 1:length(MutInfo)
            if strcmpi(MutInfo(r).Type, MutType)
                MutMotif = MutInfo(r).Motif;
                IntAfter = nt2int(MutInfo(r).NT(end));
                for k = 1:3
                    try
                        CodonLoc = codonMap(MutMotif(5-k-1:5-k+1));
                    catch
                        continue;
                    end
                    RelMutCell{CodonLoc, k}(IntAfter) = RelMutCell{CodonLoc, k}(IntAfter) + 1; 
                end
            end
        end
    end
end
MotifData.Data = RelMutCell;

varargout{1} = MotifData;
