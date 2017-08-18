%addIMGTgaps will get a V gene sequence and remove gaps. Will return a
%string data that can be used to get the gaps back.
%
%  Seq = addIMGTgaps(NoGapSeq, GapInfo)
%
%  Seq = addIMGTgaps(NoGapSeq, GapInfo, MissingLeft)
%
%  INPUT
%    NoGapSeq: nt sequence without gaps, but doesn't have to be the same
%      length as the database sequence, as long as the missing left length
%      is specfied and there are not extra 5' nts longer than the database
%      sequence.
%    GapInfo: A string sequence in the format M(-G(M(-G, where +M numbers
%      are nt counts, and -G numbers are gap counts. This MUST come from
%      the database map structure, DB.(Xmap){j, M.GapInfoLoc}. 
%      count. This is used to return a NoGapSeq to a GappedSeq
%    MissingLeft: number of nt missing 5' side of the NoGapSeq with respect
%      to the database sequence without gaps. This is used if you have an
%      incomplete V gene starting. EX: if NoGapSeq starts at the 100th nt
%      of the total V gene, then MissingLeft should be 99.
%
%  OUTPUT
%    Seq: nt sequence with IMGT "." gaps inserted according to GapInfo
%
%  EXAMPLE
%    NoGapSeq = 'caggttactgaatctgcactcac';
%    GapInfo = '3#-3#3#-3#3#-3#9#-3#5';
%    Seq = addIMGTgaps(NoGapSeq, GapInfo, 0);
%    Seq = 
%          cag...gtt...act...gaatctgca...ctcac
%
%  See also removeIMGTgaps

function Seq = addIMGTgaps(NoGapSeq, GapInfo, varargin)
%Determine and add missing nts (N) to the NoGapSeq to ensure it's the same
%length as the germline sequence.
MissingLeft = 0;
if ~isempty(varargin)
    MissingLeft = varargin{1};
end
if MissingLeft < 0
    error('%s: MissingLeft cannot be a negative number. Trim the NoGapSeq so that the 1st nt is the 1st nt of the database sequence.', mfilename)
elseif MissingLeft > 0
    NoGapSeq = [repmat('N', 1, MissingLeft) NoGapSeq];
end

%Parse the GapInfo
GapInfo = strrep(GapInfo, ' ', '');
GapMat = regexp(GapInfo, '[^\d^\-]', 'split');
for j = 1:length(GapMat)
    try
        GapMat{j} = eval(GapMat{j});
    catch
        error('%s: GapInfo contains an invalid character: %s.', mfilename, GapMat{j})
    end
end
GapMat = cell2mat(GapMat);

%Generate the final sequence with gaps
SeqCell = cell(1, length(GapMat) + 1);
s1 = 1;
for j = 1:length(GapMat)
    if GapMat(j) < 0
        SeqCell{j} = repmat('.', 1, abs(GapMat(j)));
    else
        s2 = s1 + GapMat(j) - 1;
        if s2 > length(NoGapSeq)
            s2 = length(NoGapSeq);
        end
        SeqCell{j} = NoGapSeq(s1:s2);
        s1 = s2 + 1;
    end
end
if s1 <= length(NoGapSeq)
    SeqCell{j+1} = NoGapSeq(s1:end);
end
Seq = sprintf(repmat('%s', 1, length(SeqCell)), SeqCell{:});

