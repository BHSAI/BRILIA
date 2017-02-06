%checkInputQuality will evaluate only the sequences in VDJdata to ensure it
%meets the minimum quality for BRILIA processing. Returns a logical index
%array called BadIdx that can be used to delete the bad sequences, defined
%as one with > 10% of sequence with non-nucleotide letters.
%
%  [VDJdata, BadIdx] = checkInputQuality(VDJdata,VDJheader);
%
%  See also BRILIA

function [VDJdata, varargout] = fixInputSeq(VDJdata,VDJheader) 
H = getHeaderVar(VDJheader);

BadIdx = zeros(size(VDJdata,1),1,'logical');
for j = 1:size(VDJdata,1)
    Seq = VDJdata{j,H.SeqLoc};
    
    %Make sure it is all caps
    Seq = upper(Seq);

    %Make sure no absurd characters - replace with X.
    NonNucIdx = regexp(Seq,'[^ACGTU]');
    Seq(NonNucIdx) = 'X';
    
    %Make sure NonNucIdx count is not too large, otherwise abandon Seq.
    if length(NonNucIdx)/length(Seq) > 0.10
        BadIdx(j) = 1;
        continue
    else
        VDJdata{j,H.SeqLoc} = Seq;
    end
end

if nargout == 2
    varargout{1} = BadIdx;
end
