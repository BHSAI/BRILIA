%labelRegionMEX will get a 1xN or Mx1 vector and label form 1 to N the 
%number of consecutive-matched numbers.
%    
%  Label = labelRegionMEX(Region)
%
%  INPUT
%    Region: 1xN or Mx1 vector like [0 0 1 1 2 2 2 2 0 0 1 1]
%    
%  OUTPUT
%    Label: vector of size(Region) of regions 1 to N
%
%  EXAMPLE
%    Region = [0 0 1 1 2 2 2 2 0 0 1 1];
%    Label  = labelRegionMEX(Region)
%    Label  = 
%             [1 1 2 2 3 3 3 3 4 4 5 5]
%
%
