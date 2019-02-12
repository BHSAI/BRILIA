%convProp2Charge takes 2 amino acid property strings and then returns the 
%net change in charge. Mainly looking for B(+1), A(-1), and all others (0).
%
%[ChargeArray, TotCharge] = convProp2Charge(Prop1, Prop2)
%
% INPUT
%   Prop1: starting AA property code
%   Prop2: ending AA property code
%
% OUTPUT
%   ChargeArray: 1xN double matrix of charge differences
%   TotCharge: (charges in Prop 2) - (charges in Prop 1)
%
% EXAMPLE
%   Prop1 = 'HAHHAAH' %(-3)
%   Prop2 = 'HBHHBBH' %(+3)
%   [C, D] = convProp2Charge(Prop1, Prop2)
%   C =
%       0   2   0   0   2   2   0
%   D =
%
