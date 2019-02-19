%trimMatchResultsMEX will take a binary sequence match results and trim
%matches (1's) near left, right, both, or no edges until it encounters a 
%a block where 3 out of 4 are matches (ex: 1 1 0 1). Trimmed matches are 
%0's.
% 
%  MatchResults = trimMatchResultsMEX(MatchResults, TrimSide)
%
%  INPUT
%    MatchResults: 1xN logical array
%    TrimSide ['n', 'l', 'r', 'b']: direction to trim. 'n' does nothing.
% 
%  OUTPUT
%    MatchResults: 1xN logical array with 1's turned to 0's where trimmed
%
%  EXAMPLE
%    M = [1 0 0 1 1 1 1 0 0 1] > 0;
%    TrimSide = 'b';
%    M = trimMatchResultsMEX(M, TrimSide)
%    M =
%        0   0   0   1   1   1   1   0   0   0
%
%
