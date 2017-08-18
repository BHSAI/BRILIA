%findSynMutLoc will compare a RefSeq to a Seq and identify which mutations
%are synonymous. Will also locate nonsynonymous mutations.
%
%  SynLoc = findSynMutLoc(RefSeq, Seq)0
%
%  [SynLoc, NonsynLoc] = findSynMutLoc(RefSeq, Seq)0
%
%  [SynLoc, NonsynLoc] = findSynMutLoc(RefSeq, Seq, 'Frame', Frame)
%
%  INPUT
%    RefSeq: starting NT sequence (X0)
%    Seq: resulting NT sequence (X1)
%    MutType: starting reading frame position
%
%  OUTPUT
%    SynLoc: 1xlength(Seq) logical array marking synonymous mutations
%    NonsynLoc: 1xlength(Seq) logical array marking nonsynonymous mutations
%
%  NOTE
%    If length(RefSeq) ~= length(Seq), then this will pad 0's to the
%    right side of the shorter sequence to get (Non)SynLoc with the length
%    of the longer sequence.
%
%  EXAMPLE
%    RefSeq = 'ACGTGTAGCGGG'; %TCSG
%    Seq    = 'ACCTGGCGC'; %TWRG
%    Frame  = 1;
%    [SynLoc, NonsynLoc] = findSynMutLoc(RefSeq, Seq, 'Frame', Frame);
%    SynLoc = 
%       0 0 1 0 0 0 0 0 0 0 0 0
%    NonsynLoc = 
%       0 0 0 0 0 1 1 0 0 0 0 0

function [SynLoc, varargout] = findSynMutLoc(RefSeq, Seq, varargin)
P = inputParser;
addParameter(P, 'Frame', 1, @(x) isnumerical(x))
parse(P);
Frame = P.Results.Frame;

%Ensure ambiguous letters are N
AmbigRef = regexpi(RefSeq, '[^ACGTU]');
RefSeq(AmbigRef) = 'N';
AmbigSeq = regexpi(Seq, '[^ACGTU]');
Seq(AmbigSeq) = 'N';

%Ensure RefSeq and Seq are the same lengths
if length(RefSeq) ~= length(Seq)
    AddToRef = length(Seq) - length(RefSeq);
    if AddToRef > 0
        RefSeq = sprintf('%s%s', RefSeq, Seq(length(RefSeq)+1:length(RefSeq)+AddToRef));
    else
        Seq = sprintf('%s%s', Seq, RefSeq(length(Seq)+1:length(Seq)-AddToRef));
    end
end

%Determine if there is any difference in RefSeq or Seq
SynLoc = zeros(1, length(RefSeq), 'logical');
NonsynLoc = zeros(1, length(RefSeq), 'logical');
MutLoc = RefSeq ~= Seq;
if max(MutLoc) == 0; %Nothing to mark
    if nargout == 2
        varargout{1} = NonsynLoc;
    end
    return;
end

%For each codon, check if there's different nt
for j = Frame:3:length(RefSeq)-2
    if max(MutLoc(j:j+2)) == 1 %There is a mutation, so check
        if ~isempty(strfind(RefSeq(j:j+2), 'N')); continue; end
        if ~isempty(strfind(Seq(j:j+2), 'N')); continue; end
        AA1 = nt2aa(RefSeq(j:j+2), 'ACGTonly', false);
        AA2 = nt2aa(Seq(j:j+2), 'ACGTonly', false);
        if AA1 ~= AA2 %Replacement
            NonsynLoc(j:j+2) = MutLoc(j:j+2);
        else
            SynLoc(j:j+2) = MutLoc(j:j+2);
        end
    end
end

if nargout == 2
    varargout{1} = NonsynLoc;
end
