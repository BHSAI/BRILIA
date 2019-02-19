%cmprSeqMEX will compare two sequences and return a logical array where
%1 = match and 0 = mismatch. 
%  
%  Match = cmprSeqMEX(SeqA, SeqB, Alphabet)
%   
%  INPUT
%    SeqA: 1st sequence (CASE SENSITIVE!)
%    SeqB: 2nd sequence (CASE SENSITIVE!)
%    Alphabet ['n', 'a', 'r']: nucleotide, amino acid, or random sequenc
%      'n': X and N are wildcard for DNA/RNA (default)
%      'a': X is wildcard for AA
%      'r': no wildcards for character
%    
%  OUTPUT
%    Match: 1-row logical array of matches (1) and mismatch(0)
%
%  EXAMPLE
%    SeqA = 'ACGTXXNNACGT';
%    SeqB = 'ACGTACGTACGT';
%
%    Match = cmprSeqMEX(SeqA, SeqB, 'n')
%        =  1  1  1  1  1  1  1  1  1  1  1  1
%    Match = cmprSeqMEX(SeqA, SeqB, 'a')
%        =  1  1  1  1  1  1  0  0  1  1  1  1
%    Match = cmprSeqMEX(SeqA, SeqB, 'r')
%        =  1  1  1  1  0  0  0  0  1  1  1  1
%
%
