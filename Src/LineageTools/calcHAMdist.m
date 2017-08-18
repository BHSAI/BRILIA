%calcHAMdist will do pairwise comparison of two sequences and return the
%number of mismatched letters. By default, an 'X' is a wildcard.
%
%  HamDist = calcHAMdist(Seq1,Seq2)
%
%  HamDist = calcHAMdist(Seq1,Seq2,WildCard)
%
%  INPUT
%    Seq1: sequence string
%    Seq2: Mx1 cell of sequences, or a sequence string
%    WildCard ['y' 'n']: Determines if an 'X' character is a wildcard or
%      not. Default is WildCard = 'y'.
%
%  OUTPUT
%    HamDist: Mx1 matrix of mismatched letters between Seq1 and Seq2.
%
%  NOTE
%    'X' is a wildcard match.
%
%    If Seq1 and Seq2 are variable lengths, then HamDist is calculated as
%    number of mismatched letters up to the length of the shorter sequence,
%    plus the difference in sequence lengths:
%
%  EXAMPLES
%    Seq1 = 'ACXTAT' 
%    Seq2 = 'ATGTGT'
%    HamDist = calcHAMdist(Seq1,Seq2)
%    HamDist = 
%               2
%
%    HamDist = calcHAMdist(Seq1,Seq2,'n')
%    HamDist = 
%               3
%
%    Seq2 = {'ACXT','ATGTGT'};
%    HamDist = calcHAMdist(Seq1,Seq2,'n')
%    HamDist = 
%               2
%               3
function Dist = calcHAMdist(Seq1,Seq2,varargin)
%Determine wild card X matching
MatchWild = 'y'; %Default, yes
if ~isempty(varargin)
    MatchWild = lower(varargin{1}(1));
    if ~ismember(MatchWild,{'y','n'})
        error('calcHAMdist: WildCard y or n input is invalid');
    end
end

%Need Seq1 to be char
if iscell(Seq1); 
    Seq1 = Seq1{1};
end

%Want Seq2 to be cell 
if ischar(Seq2) 
    Seq2 = {Seq2};
end

%Calc Hamming distances without wildcard X matching enabled
Dist = zeros(length(Seq2),1);
if MatchWild == 'n'
    for j = 1:length(Seq2)
        if length(Seq1) == length(Seq2{j})
            SeqMatch = (Seq1 == Seq2{j});
            Dist(j) = sum(~SeqMatch);
        else
           SeqLens = sort([length(Seq1) length(Seq2{j})]);
           SeqMatch = Seq1(1:SeqLens(1)) == Seq2{j}(1:SeqLens(1));
           Dist(j) = sum(~SeqMatch) + diff(SeqLens);        
        end
    end
    
%Calc Hamming distances with wildcard X matching enabled
elseif MatchWild == 'y'
    Xloc1 = Seq1 == 'X';
    for j = 1:length(Seq2)
        Xloc2 = Seq2{j} == 'X';
        if length(Seq1) == length(Seq2{j})
            SeqMatch = (Seq1 == Seq2{j}) | Xloc1 | Xloc2;
            Dist(j) = sum(~SeqMatch);
        else
            SeqLens = sort([length(Seq1) length(Seq2{j})]);
            SeqMatch = Seq1(1:SeqLens(1)) == Seq2{j}(1:SeqLens(1)) | Xloc1(1:SeqLens(1)) | Xloc2(1:SeqLens(1)); 
            Dist(j) = sum(~SeqMatch) + diff(SeqLens);        
        end
    end
end
