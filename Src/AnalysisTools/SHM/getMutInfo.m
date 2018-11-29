%getMutInfo will compare two DNA/RNA sequences and return the details of
%the mutations, including X -> Y, Nth position, nnXnn 5-letter motifs
%around the initial sequence's mutant, and synonymous mutations.
%
%  MutInfo = getMutInfo(SeqA, SeqB)
%
%  MutInfo = getMutInfo(SeqA, SeqB, Frame)
%
%  MutInfo = getMutInfo(..., 'cell')
%
%  INPUT
%    SeqA: starting DNA/RNA sequence
%    SeqB: resulting DNA/RNA sequence
%    Frame [1-3]: reading frame number (required for synonymous mutation).
%      If empty, .AA and .IsSyn fields will also have empty values.
%    'cell': returns output as a structure or cell format
%
%  OUTPUT (default)
%    MutInfo: nonscalar structure of each MUTATED nt
%      .Idx    - the nth position of nt mutation
%      .NT     - char like 'x>y' indicating nucleotide x to y mutation
%      .AA     - char like 'X>Y' indicating amino acid X to Y mutation
%      .Motif  - 5-letter motif around mutant (X = unknown/filler) nnXnn
%      .IsSyn  - 1 or 0 for synonymous mutations (requires Frame input)
%
%  OUTPUT ('cell' option)
%    MutInfo: scalar structure
%      .Header - 1xN cell of data header strings
%      .Data   - MxN cell of data
%
%  NOTE
%    If adjacent nts are mutated AND lead to a nonsynon. mutation,
%    then both nts are marked nonsynon. even though a single change will
%    lead to a syn. mutation.
%
%  EXAMPLE
%    SeqA =  'ACGTGATGATACGTGTACCCATGAA'; %RDDTCTHE,    SynLoc 0000100000000000000010000
%    SeqB =  'ATGTAATGACACGTGTGCCCGCGTA'; %CNDTCARV, NonsynLoc 0100000001000000100001010
%    Frame = 2;
%    MutInfo = getMutInfo(SeqA, SeqB, Frame)
%    MutInfo(1) = 
%         Idx: 2
%          NT: 'C > T'
%          AA: 'R > C'
%       Motif: 'NACGT'
%       IsSyn: 1
% 
%    MutInfo = getMutInfo(SeqA, SeqB, Frame, 'cell')
%    MutInfo.Header = 
%       'Idx'   'NT'       'AA'       'Motif'    'IsSyn'
%    MutInfo.Data = 
%       [ 2]    'C > T'    'R > C'    'NACGT'    [1]
%       [ 5]    'G > A'    'D > N'    'GTGAT'    [1]
%       [10]    'T > C'    'D = D'    'GATAC'    [0]
%       [17]    'A > G'    'T > A'    'GTACC'    [1]
%       [21]    'A > G'    'H > R'    'CCATG'    [1]
%       [22]    'T > C'    'H > R'    'CATGA'    [1]
%       [24]    'A > T'    'E > V'    'TGAAN'    [1]

function Out = getMutInfo(SeqA, SeqB, varargin)
%Check for 'cell' format option
UseCellFmt = false;
for j = 1:length(varargin)
    if ischar(varargin{j}) && strcmpi(varargin{j}, 'cell')
        UseCellFmt = true;
        varargin(j) = [];
        break
    end
end

%Check for Frame input
Frame = 0;
if ~isempty(varargin) 
    if isnumeric(varargin{1}) && ~isempty(intersect(varargin{1}, [1:3]))
        Frame = varargin{1};
    else
        error('%s: Check the Frame input.', mfilename);
    end
end

%Get the nucleotide mutation info
MutIdx = find(~cmprSeqMEX(SeqA, SeqB, 'n'));
MutNum = length(MutIdx);
MutInfo(1:MutNum) = struct('Idx', 0, 'NT', '', 'AA', '', 'Motif', 'NNNNN', 'IsSyn', []);
for j = 1:length(MutIdx)
    MutInfo(j).Idx = MutIdx(j);
    MutInfo(j).NT = [SeqA(MutIdx(j)) ' > ' SeqB(MutIdx(j))];
    
    %Determine Motif
    L = MutIdx(j) - 2;
    R = MutIdx(j) + 2;
    AddL = 0;
    AddR = 0;
    if L < 1
        AddL = 1 - L;
    end
    if R > length(SeqA)
        AddR = R - length(SeqA);
    end
    MutInfo(j).Motif(1+AddL:5-AddR) = SeqA(L+AddL:R-AddR);
    
    %Determine Syn/Nonsyn and AA1>AA2
    if Frame > 0
        CodonNum = ceil((MutIdx(j) - (Frame-1)) / 3);
        S = (Frame-1) + 1 + 3*(CodonNum-1);
        E = (Frame-1) + 3*CodonNum;
        
        if S < 1 || E > length(SeqA); continue; end
        NT = [SeqA(S:E) SeqB(S:E)];
        AA = nt2aa(regexprep(NT, 'X|x', 'N'), 'acgtonly', false, 'alternative', false);
        if AA(1) == AA(2)
            MutInfo(j).AA = [AA(1) ' = ' AA(2)];
            MutInfo(j).IsSyn = 0;
        else
            MutInfo(j).AA = [AA(1) ' > ' AA(2)];
            MutInfo(j).IsSyn = 1;
        end
    end
end

%Convert to cell format if needed
if UseCellFmt
    Out.Header = fieldnames(MutInfo)';
    Out.Data = squeeze(struct2cell(MutInfo))';
else
    Out = MutInfo;
end