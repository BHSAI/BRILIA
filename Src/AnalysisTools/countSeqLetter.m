%countSeqLetter will take a vector of NT or AA sequences, and then count
%how many letters are there per position. Returns either a 20xN amino acid
%array (21st is ambiguous or special), or a 4xN nt array. Unique characters
%are ignored by default.
%
%  [Freq, Letter] = countSeqletter(Seq)
%
%  [Freq, Letter] = countSeqletter(Seq, Alphabet)
%
%  INPUT
%    Seq: Nx1 cell array or NxM char array of sequences.
%    Alphabet: 'NT' or 'AA' specifies which type of sequences to be used.
%   
%  OUTPUT
%    Freq: 4xM or 20xM matrix counting the frequence of each nucleotide.
%      Ignores ambiguous characters.
%    Letter: 4x1 or 20x1 cell listing the order of the NT or AA letter
%      returned by countSeqLetter.
%
%  EXAMPLE
%    Seq = {'ACGUGACGTG';
%           'CCGUGAAGTT';
%           'TCCUGAAGTT'};
%    [Freq, Letter] = countSeqLetter(Seq)
%     Freq =
%        1   0   0   0   0   3   2   0   0   0
%        1   3   1   0   0   0   1   0   0   0
%        0   0   2   0   3   0   0   3   0   1
%        1   0   0   3   0   0   0   0   3   2
%     Letter = 
%       'A'
%       'C'
%       'G'
%       'T'
function varargout = countSeqLetter(Seq, varargin)
%Convert to character first
if iscell(Seq)
    Seq = char(Seq);
end

%Determine if it's an NT or AA
if isempty(varargin) %Autodetect
    NonNTLoc = regexpi(Seq(:)','[^ACGTU]','once');
    if ~isempty(NonNTLoc) 
        Alphabet = 'AA';
    else
        Alphabet = 'NT';
    end
elseif ischar(varargin{1})
    Alphabet = upper(varargin{1});
    if isempty(regexp(Alphabet,'NT|AA','once'))
        Alphabet = 'AA';
        warning('countSeqLetter: Invalid Alphabet name. Using ''AA''');
    end
end

%Count how many letters there are
if strcmpi(Alphabet,'AA')
    Counts = zeros(20,size(Seq,2));
    for j = 1:size(Seq,2)
        Stat = aacount(Seq(:,j),'Ambiguous','Ignore');
        Counts(:,j) = cell2mat(struct2cell(Stat));
    end
else
    Counts = zeros(4,size(Seq,2));
    for j = 1:size(Seq,2)
        Stat = basecount(Seq(:,j),'Ambiguous','Ignore');
        Counts(:,j) = cell2mat(struct2cell(Stat));
    end
end

%Return counts and letters
if nargout >= 1
    varargout{1} = Counts;
    if nargout >= 2
        varargout{2} = fieldnames(Stat);
    end
end
