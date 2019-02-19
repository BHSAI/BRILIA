%alignSeqMEX will align sequence SeqA and SeqB using the BRILIA alignment
%protocol. Consecutive matches add C^2 score and consecutive mismatches
%subtract M^2 score, where C and M is the length of each contiguous
%segment. If AllowedMiss is specified, then C segments can be elongated
%only if there is a 1-letter miss, and if the number of AllowedMiss has not
%been reached, going left to right. More information is provided in the
%reference for BRILIA.
%
%  Score = alignSeqMEX(SeqA, SeqB)
%
%  [Score, StartAt] = alignSeqMEX(SeqA, SeqB)
%
%  [Score, StartAt, MatchAt] = alignSeqMEX(SeqA, SeqB)
%
%  [Score, StartAt, MatchAt, Alignment] = alignSeqMEX(SeqA, SeqB)
%
%  ... = alignSeqMEX(SeqA, SeqB, MissRate, Alphabet, ExactMatch, TrimSide,
%          PenaltySide, PreferSide)
%  
%  INPUTS
%    SeqA: Character sequence. X = wildcard match. Z = do not match.
%    SeqB: Character sequence. X = wildcard match. Z = do not match.
%    MissRate [Integer >= 0]: How many missed nt/bp are allowed for
%      elongating a contiguous C region for alignment score. Ex: MissRate
%      0.02 = 2 of shorter seq length.
%    Alphabet ['n' 'a' 'r']: Specify nt, aa, or random char. For nt, N and
%      X char are wildcard. for aa, only X is wildcard. For random, no
%      wildcards.
%    ExactMatch ['n' 'y']: If SeqA and SeqB are same lengths and you
%      want to do a simple, direct alignment, then use 'y'.
%    TrimSide ['n' 'l','r' 'b']: Will perform a 3/4 match on the side:
%      none, left, right, or both.
%    PenaltySide ['n' 'l' 'r']: Will reduce alignment score if
%      none, left, or right portion of SeqA overhangs SeqB and is unused.
%      Penalty is (unused nt)^2.
%    PreferSide ['n' 'l' 'r']: If a tie alignment is found, will favor the
%      alignment closers to this side of SeqA.
%
%  OUTPUTS
%    Score: 2x1 matrix, Score(1) is the # of hits, and Score(2) is the
%      aligment score of BRILIA.
%    StartAt: 2x1 matrix, StartAt(1) is always 1 the SeqA start, and
%      StartAt(2) is the SeqB start relative to the nth letter of SeqA. 
%      EX: if StartAt(2) = 3, then SeqB aligns with A if its 1st letter
%      start below the 3rd letter of SeqB. If StartAt(2) = -3, then the 1st
%      letter of SeqB start left of SeqA's 1st letter by 3 letters (missing
%      3 letters).
%    MatchAt: 2x1 matrix showing the 1st and last location of the '|' in
%      the Alignment char matrix.
%    Alignment: 3xM char matrix showing the alignment results of SeqA and
%      SeqB. '|' in the 2nd char row marks locations of matches, while ' '
%      in the 1st and 3rd row are unmatched letters.
%
%--------------------------------------------------------------------------
%  EXAMPLES
%    Case1) "TrimSide" edge cleaning
%      SeqA = 'TAATAATTAAT'
%      SeqB = 'TCCTAATTGGT'
%      [Score, StartAt, MatchAt, Alignment] = alignSeqMEX(SeqA, SeqB)
%      Score =
%            7
%           19
%      Alignment =
%         TAATAATTAAT
%         |  |||||  |
%         TCCTAATTGGT
%      StartAt =
%            1
%            1
%      MatchAt =
%            1
%           11
%
%      [Score, ~, ~, Alignment] = alignSeqMEX(SeqA, SeqB, 0, 'n', 'n', 'l')
%       Alignment =
%         TAATAATTAAT
%            |||||  |
%         TCCTAATTGGT
%
%      [Score, ~, ~, Alignment] = alignSeqMEX(SeqA, SeqB, 0, 'n', 'n', 'r')
%       Alignment =
%         TAATAATTAAT
%         |  |||||   
%         TCCTAATTGGT
%
%      [Score, ~, ~, Alignment] = alignSeqMEX(SeqA, SeqB, 0, 'n', 'n', 'b')
%       Alignment =
%         TAATAATTAAT
%            |||||   
%         TCCTAATTGGT
%
%    Case2) "PreferSide" tie-breaking
%      SeqA = 'AAAA'
%      SeqB = 'CAAAATTAAAATTAAAAC'
%      [~, ~, ~, Alignment] = alignSeqMEX(SeqA, SeqB, 0, 'n', 'n', 'n', 'n', 'l')
%       Alignment =
%         -AAAA-------------
%          ||||             
%         CAAAATTAAAATTAAAAC
%
%      [~, ~, ~, Alignment] = alignSeqMEX(SeqA, SeqB, 0, 'n', 'n', 'n', 'n', 'r')
%       Alignment =
%         -------------AAAA-
%                      |||| 
%         CAAAATTAAAATTAAAAC
%
%      [~, ~, ~, Alignment] = alignSeqMEX(SeqA, SeqB, 0, 'n', 'n', 'n', 'n', 'n')
%       Alignment =
%         -------AAAA-------
%                ||||       
%         CAAAATTAAAATTAAAAC
%
%    Case3) Wildcard matching and exact matching
%      SeqA = 'CGAANCAA'
%      SeqB = 'ACGAACGA'
%      [~, ~, ~, Alignment] = alignSeqMEX(SeqA, SeqB, 0, 'n', 'n')
%       Alignment =
%         -CGAAXCAA
%          ||||| |  
%         ACGAACGA-
%
%      [~, ~, ~, Alignment] = alignSeqMEX(SeqA, SeqB, 0, 'n', 'y')
%       Alignment =
%         CGAAXCAA
%            ||| |
%         ACGAACGA
%
%    Case4) AllowedMiss scoring changes
%      SeqA = 'ACGTGGTA'
%      SeqB = 'ACATGATA'
%      [Score, ~, ~, Alignment] = alignSeqMEX(SeqA, SeqB, 0, 'n', 'y')
%       Score =
%           6
%          10
%       Alignment =
%         ACGTGGTA
%         || || ||
%         ACATGATA
%
%      [Score, ~, ~, Alignment] = alignSeqMEX(SeqA, SeqB, 1/8, 'n', 'y')
%        Score =
%            6
%           19
%
%
%
