%fixInputSeq will check all input sequences within VDJdata to ensure it
%meets the minimum quality for BRILIA processing. Returns a logical index
%array called BadIdx that can be used to delete the bad sequences, which
%are those with > 10% of sequence with non-nucleotide letters.
%
%  VDJdata = fixInputSeq(VDJdata, Map)
%
%  [VDJdata, BadIdx] = fixInputSeq(VDJdata, Map, MaxErrPerc)
%
%  INPUT
%    VDJdata: main BRILIA data cell
%    Map: structure of BRILIA data header index
%    MaxErrRate: max fraction of nts in seq that is tolerable (0.1 default)
%
%  OUTPUT
%    VDJdata: modified VDJdata with corrected sequences
%    BadLoc: Nx1 logical array of entries have > 10% non-nt letters

function [VDJdata, varargout] = fixInputSeq(VDJdata, Map, varargin)
SeqIdx = [Map.hSeq; Map.lSeq];
SeqIdx(SeqIdx == 0) = [];

%Set the seq error %, which is the # of non-nt letters / length of seq
MaxErrRate = 0.1;
if ~isempty(varargin) && isnumeric(varargin)
    MaxErrRate = varargin{1};
end

%Replace ambig seq nts with X, and mark those with > MaxErrPerc in BadIdx
BadLoc = zeros(size(VDJdata, 1), 1, 'logical');
for j = 1:size(VDJdata, 1)
    for k = 1:length(SeqIdx)
        %Make uppercase, and replace absurd characters with X.
        VDJdata{j, SeqIdx(k)} = upper(VDJdata{j, SeqIdx(k)});
        NonNucIdx = regexp(VDJdata{j, SeqIdx(k)}, '[^ACGTU]');
        VDJdata{j, SeqIdx(k)}(NonNucIdx) = 'X';

        %Make sure NonNucIdx count is not too large, otherwise abandon Seq.
        if length(NonNucIdx)/length(VDJdata{j, SeqIdx(k)}) > MaxErrRate
            BadLoc(j) = 1;
        end
    end
end

if nargout == 2
    varargout{1} = BadLoc;
end