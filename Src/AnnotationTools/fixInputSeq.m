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
if nargin < 3
    MaxErrRate = 0.1;
end
SeqIdx = [Map.hSeq; Map.lSeq];
SeqIdx = SeqIdx(SeqIdx > 0);
VDJdata(:, SeqIdx) = regexprep(upper(VDJdata(:, SeqIdx)), '[^ACGTU]', 'N');
MissRate = cellfun(@(x) sum(x == 'N') / length(x), VDJdata(:, SeqIdx));
BadLoc = any(MissRate > MaxErrRate, 2);

FunctIdx = [Map.hFunct Map.lFunct];
FunctIdx = FunctIdx(FunctIdx > 0);
VDJdata(BadLoc, FunctIdx) = {'I'};