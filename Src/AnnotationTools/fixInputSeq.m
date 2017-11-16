%fixInputSeq will check all input sequences within VDJdata to ensure it
%meets the minimum quality for BRILIA processing. Returns a logical index
%array called BadIdx that can be used to delete the bad sequences, which
%are those with > 10% of sequence with non-nucleotide letters.
%
%  VDJdata = fixInputSeq(VDJdata, VDJheader)
%
%  [VDJdata, BadIdx] = fixInputSeq(VDJdata, VDJheader)
%
%  [VDJdata, BadIdx] = fixInputSeq(VDJdata, VDJheader, MaxErrPerc)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    VDJheader: main BRILIA header cell
%    MaxErrRate: max fraction of nts in seq that is tolerable (0.1 default)
%
%  OUTPUT
%    VDJdata: modified VDJdata with corrected sequences
%    BadIdx: Nx1 logical array of entries have > 10% non-nt letters

function [VDJdata, varargout] = fixInputSeq(VDJdata, Map, varargin)
%Extract all sequence locations
SeqLoc = [Map.hSeq; Map.lSeq];
SeqLoc(SeqLoc == 0) = [];

%Determine max error %
MaxErrRate = 0.1;
if ~isempty(varargin) && isnumeric(varargin)
    MaxErrRate = varargin{1}(1);
end

%Replace ambig seq nts with X, and mark those with > MaxErrPerc in BadIdx
BadIdx = zeros(size(VDJdata, 1), 1, 'logical');
for j = 1:size(VDJdata, 1)
    for k = 1:length(SeqLoc)
        %Make uppercase, and replace absurd characters with X.
        Seq = upper(VDJdata{j, SeqLoc(k)});
        NonNucIdx = regexp(Seq, '[^ACGTU]');
        Seq(NonNucIdx) = 'X';

        %Make sure NonNucIdx count is not too large, otherwise abandon Seq.
        if length(NonNucIdx)/length(Seq) > MaxErrRate
            BadIdx(j) = 1;
        else
            VDJdata{j, SeqLoc(k)} = Seq;
        end
    end
end

if nargout == 2
    varargout{1} = BadIdx;
end
