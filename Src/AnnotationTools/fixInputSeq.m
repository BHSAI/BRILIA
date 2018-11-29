%fixInputSeq will check all input sequences within VDJdata to ensure it
%meets the minimum quality for BRILIA processing. Returns a logical index
%array called BadIdx that can be used to delete the bad sequences, which
%are those with > 10% of sequence with wildcard or non-ACGTU letters.
%
%  [VDJdata, BadLoc] = fixInputSeq(VDJdata, Map)
%
%  [VDJdata, BadLoc] = fixInputSeq(VDJdata, Map, MaxErrPerc)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    Map: structure of BRILIA data header index
%    MaxErrRate: max fraction of nts in seq that is tolerable (0.1 default)
%
%  OUTPUT
%    VDJdata: modified VDJdata with corrected sequences
%    BadLoc: Nx1 logical array of entries have > 10% non-nt letters
%
function [VDJdata, BadLoc] = fixInputSeq(VDJdata, Map, MaxErrRate)
BadLoc = repelem(false, size(VDJdata,1), 1);
if isempty(VDJdata); return; end
if nargin < 3
    MaxErrRate = 0.1;
end

SeqIdx = nonzeros([Map.hSeq; Map.lSeq]);
VDJdata(:, SeqIdx) = regexprep(upper(VDJdata(:, SeqIdx)), '[^ACGTU]', 'N');

BadLoc = repelem(false, size(VDJdata, 1), 1);
for c = 1:length(Map.Chain)
    C = lower(Map.Chain(c));
    SeqIdx = Map.([C 'Seq']);
    OverSeq5 = Map.([C 'OverSeq5']);
    OverSeq3 = Map.([C 'OverSeq3']);
    for j = 1:size(VDJdata, 1)
        S = find(VDJdata{j, SeqIdx} ~= 'N', 1, 'first');
        E = find(VDJdata{j, SeqIdx} ~= 'N', 1, 'last');
        VDJdata{j, OverSeq5} = VDJdata{j, SeqIdx}(1:S-1);
        VDJdata{j, OverSeq3} = VDJdata{j, SeqIdx}(E+1:end);
        VDJdata{j, SeqIdx} = VDJdata{j, SeqIdx}(S:E);
        MissRate = sum(VDJdata{j, SeqIdx} == 'N') / length(VDJdata{j, SeqIdx});
        BadLoc(j) = BadLoc(j) | (MissRate > MaxErrRate);
    end
end

FunctIdx = nonzeros([Map.hFunct Map.lFunct]);
VDJdata(BadLoc, FunctIdx) = {'I'};