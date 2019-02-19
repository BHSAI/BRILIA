%countHotspotsMEX will count the number of hotspot motifs in a sequence. 
%Partial motifs at the sequence edge will NOT be considered a hotspot.
%
%  N = countHotspotsMEX(Seq)
% 
%  INPUT
%    Seq: nt sequence string or cell of sequence strings (all CAPS)
%
%  OUTPUT
%    N: scalar or matrix of number of hotspots
%
%  EXAMPLE
%
%    Seq = {'TGCTGCTGC';
%           'TAGGUAG'};
%    [N, S] = countHotspotsMEX(Seq);
%
%
%
