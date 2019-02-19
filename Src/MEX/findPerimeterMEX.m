%findPerimeterMEX will get the row, column, and linear index of pixels 
%around a certain pixel. This is faster than using MATLAB if/else statements
%and sub2ind functions.
%
%  [R, C, Idx] = findPerimeterMEX(Size, Rc, Cc)
%
%  INPUT
%    Size: size of the 2D matrix
%    Rc: row number of the pixel
%    Cc: col number of the pixel   
% 
%  OUTPUT
%    R: row index of the perimeter pixels
%    C: column index of the perimeter pixels
%    Idx: linear index of the perimeter pixels
%
%  EXAMPLE
%    G = rand(10);
%    [R, C, Idx] = findPerimeterMEX(size(G), 1, 1)
%    R = 
%        2
%        1
%        2
%    C = 
%        1
%        2
%        2
%    Idx =
%        2
%       11
%       12
%
%
%
