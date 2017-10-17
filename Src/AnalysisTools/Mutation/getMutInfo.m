%getMutInfo will compare two sequences and return the details of the
%mutations, including X -> Y, position, +2 X -2 letter motifs of the
%RefSeq, S synonymous or N nonsynonymous mutation. This is the central
%function for MutationAnalysis.
%
%  MutInfo = getMutInfo(RefSeq, Seq)  
%
%  MutInfo = getMutInfo(RefSeq, Seq, 'Frame', Frame, ...)  
%
%  MutInfo = getMutInfo(RefSeq, Seq, 'ReturnAs', ReturnAs, ...)  
%
%  INPUT
%    RefSeq: reference or ancestor sequence
%    Seq: descendant sequence
%    Frame: reading frame number (required for S / N determination). If
%      empty, does not determine syn/nonsyn mutations.
%    ReturnAs: 'struct' or 'cellstruct' output format
%      'struct': MutInfo have fields Before, After, Motif, Type, Idx
%      'cellstruct': MutInfo have fields Header, Data
%
%  OUTPUT
%    
%
%  WARNING!
%    RefSeq and Seq should be same lower or upper case, otherwise it will
%    be mismatched. 'a not equal to 'A'
%
%  EXAMPLE
%    Frame = 2
%    RefSeq =  'ACGTGATGATACGTGTACCCATGAA' %RDDTCTHE
%    Seq =     'ATGTAATGACACGTGTGCCCGCGTA' %CNDTCARV
%   %SynLoc     0000100000000000000010000
%   %NonsynLoc  0100000001000000100001010
%    MutInfo = getMutInfo(RefSeq, Seq, 'Frame', Frame)
%    struct2cell(MutInfo)
%    ans(:,:,1) = 
%             'C'
%             'T'
%             'XACGT'
%             'N'
%             [    2]
% 
%    ans(:,:,2) = 
%             'G'
%             'A'
%             'GTGAT'
%             'S'
%             [    5]
%
%    MutInfo = getMutInfo(RefSeq, Seq, 'Frame', Frame, 'ReturnAs', 'cellstruct')
%    MutInfo.Header = 
%      'Before' 'After' 'Motif'   'Type'  'Idx'
%    MutInfo.Data = 
%         'C'    'T'    'XACGT'    'N'    [ 2]
%         'G'    'A'    'GTGAT'    'S'    [ 5]
%         'T'    'C'    'GATAC'    'N'    [10]
%         'A'    'G'    'GTACC'    'N'    [17]
%         'A'    'G'    'CCATG'    'S'    [21]
%         'T'    'C'    'CATGA'    'N'    [22]
%         'A'    'T'    'TGAAX'    'N'    [24]
%    
function MutInfo = getMutInfo(RefSeq, Seq, varargin)
Frame    = parseInputVar('Frame', 0, @(x) isnumeric(x) && x >= 0 && x <= length(RefSeq), varargin{:});
ReturnAs = parseInputVar('ReturnAs', 'struct', @(x) any(strcmpi(x, {'struct', 'cellstruct'})), varargin{:});

%Ensure everything is upper case
RefSeq = upper(RefSeq);
RefSeq(regexp(RefSeq, '[^ACGTU]')) = 'N';
Seq = upper(Seq);
Seq(regexp(RefSeq, '[^ACGTU]')) = 'N';

%Ensure RefSeq and Seq are the same lengths
if length(RefSeq) ~= length(Seq)
    AddToRef = length(Seq) - length(RefSeq);
    if AddToRef > 0
        RefSeq = sprintf('%s%s', RefSeq, Seq(length(RefSeq)+1:length(RefSeq)+AddToRef));
    else
        Seq = sprintf('%s%s', Seq, RefSeq(length(Seq)+1:length(Seq)-AddToRef));
    end
end

if Frame > 0
    [SynLoc, NonsynLoc] = findSynMutLoc(RefSeq, Seq, 'Frame', Frame);
    MutIdx = find(SynLoc | NonsynLoc);
else
    MutIdx = find(~(RefSeq == Seq | Seq == 'X' | RefSeq == 'X'));
end

%Get the mutation details
MutNum = length(MutIdx);
MutInfo(1:MutNum) = struct('Before', '', 'After', '', 'Motif', '', 'Type', '', 'Idx', 0);
for j = 1:length(MutIdx)
    %Fill in basic data
    MutInfo(j).Idx = MutIdx(j);
    MutInfo(j).Before = RefSeq(MutIdx(j));
    MutInfo(j).After = Seq(MutIdx(j));
    
    %Determine Motif
    s1 = MutIdx(j) - 2;
    AddLeftMotif = '';
    if s1 < 1
        AddLeftMotif = repmat('X', 1, abs(s1 - 1));
        s1 = 1;
    end
    s2 = MutIdx(j) + 2;
    AddRightMotif = '';
    if s2 > length(RefSeq)
        AddRightMotif = repmat('X', 1, s2 - length(RefSeq));
        s2 = length(RefSeq);
    end
    MutInfo(j).Motif = cat(2, AddLeftMotif, RefSeq(s1:s2), AddRightMotif);
    
    %Determine if it's silent or not
    if Frame > 0
        if SynLoc(MutIdx(j))
            MutInfo(j).Type = 'S';
        else
            MutInfo(j).Type = 'N';
        end
    end
end

%See if you should convert to structure with format A.header, A.data
if strcmpi(ReturnAs, 'cellstruct')
    MutHeader = fieldnames(MutInfo)';
    MutData = cell(length(MutInfo), length(MutHeader));
    for j = 1:length(MutInfo)
        MutData(j, :) = struct2cell(MutInfo(j))';
    end
    clear MutInfo;
    MutInfo.Header = MutHeader;
    MutInfo.Data = MutData;
end
