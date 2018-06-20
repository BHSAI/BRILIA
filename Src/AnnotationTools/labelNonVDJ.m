%removeNonVDJ will look evaluate if sequences are too short or have too
%many errors in the V / J gene to be considered a valid sequence. This is
%used to filter out total RNAseq data that sequenes non-VDJ segements.
%
%  [VDJdata, DelThis] = labelNonVDJ(VDJdata, Map)
%
%  INPUT
%    VDJdata: BRILIA data table
%    Map: VDJheader or structure map of data column
%    ThreshHold [0 to 1]: Hamming distance fraction threshold for
%      considering a sequence is bad. 0.4 would mean if 40% of sequence in 
%      V and J are mismatched with germline, then delete. HamDiff >= 0.40,
%      then mark as invalid I.
%
%  OUTPUT
%    VDJdata: Same as input, but the H-Functional and L-Functional have "I"
%      for sequences in which the 
%    DelThis: logical index indicating where the changes have been made

function [VDJdata, DelLoc] = labelNonVDJ(VDJdata, Map, ThreshHold)

error('%s is obsolete. Use labelSeqQuality instead', mfilename);


if nargin < 3
    ThreshHold = 0.4;
end

if iscell(Map)
    Map = getVDJmapper(Map);
elseif ~isstruct(Map)
    error('%s: The 2nd input, Map, should be a struct or cell storing header info.', mfilename);
end

DelLoc = zeros(size(VDJdata, 1), 1, 'logical');
for c = 1:length(Map.Chain)
    Chain = lower(Map.Chain(c));
    SegLen = VDJdata(:, Map.([Chain 'Length']));
    Seq = VDJdata(:, Map.([Chain 'Seq']));
    RefSeq = VDJdata(:, Map.([Chain 'RefSeq']));
    
    for j = 1:length(Seq)
        if DelLoc(j); continue; end
        if isempty(SegLen{1}) || isempty(SegLen{end}); continue; end
        if isempty(Seq{j}) || isempty(RefSeq{j}); continue; end
        CmprLoc = zeros(1, length(Seq{j}), 'logical');
        CmprLoc(1:SegLen{j, 1}) = 1;
        CmprLoc(end-SegLen{j, end}+1:end) = 1;
        Score = alignSeqMEX(Seq{j}(CmprLoc), RefSeq{j}(CmprLoc), 0, 'n', 'y', 'n', 'n', 'n');
        if Score(1)/length(RefSeq{j}) < (1 - ThreshHold)
            DelLoc(j) = 1;
            continue
        end
    end
end

VDJdata(DelLoc, Map.([Chain 'Funct'])) = {'I'}; %Invalid
if any(DelLoc)
    fprintf('  %s: %d invalid sequence(s) marked.\n', mfilename, sum(DelLoc));
end