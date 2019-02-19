%findConsecOnesMEX will get a 1xN or Mx1 vector and label form 1 to N the
%number of consecutive-matched ones.
%
%  Label = findConsecOnesMEX(Region)
%    
%  INPUT
%    Region: 1xN or Mx1 vector like [0 0 1 1 1 1 0 0 0 0 1 1]
%    
%  OUTPUT
%    Label: vector of size(Region) of regions 1 to N
%
%  EXAMPLE
%    Region = [0 0 1 1 1 1 0 0 0 0 1 1];
%    Label  = findConsecOnesMEX(Region)
%    Label  = 
%             [0 0 1 1 1 1 0 0 0 0 2 2]
%
%
