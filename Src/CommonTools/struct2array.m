%struct2array description
%
%  C = struct2array(S, RowOrCol)
%
%  INPUT
%    S: nonscalor structure
%    RowOrCol ['row' 'col' 'input']: specify if the put each row as a field and 
%      data, or each col as a field and data. Default is 'col'. 
%
%  EXAMPLE
%    S(1:3) = struct('a', 1, 'b', 2, 'c', 'txt');
%    C = struct2array(S, 'row')
%    C = 
%        'a'    [  1]    [  1]    [  1]
%        'b'    [  2]    [  2]    [  2]
%        'c'    'txt'    'txt'    'txt'
%
%    C = struct2array(S, 'col')
%        'a'    'b'    'c'  
%        [1]    [2]    'txt'
%        [1]    [2]    'txt'
%        [1]    [2]    'txt'

function C = struct2array(S, RowOrCol)
F = fieldnames(S);
C = squeeze(struct2cell(S));
if nargin < 2
    RowOrCol = 'col';
end
switch lower(RowOrCol)
    case 'row'
        C = horzcat(F, C);
    case 'col'
        C = vertcat(F', C');
    case 'input'
        C = vertcat(F', C');
        C = C(:)';
    otherwise
        error('%s: Unrecognized option "%s". Use ''row'' or ''col''.', mfilename, RowOrCol);
end