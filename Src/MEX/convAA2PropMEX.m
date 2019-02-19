%convAA2PropMEX will translate AA seq to its property code. Input is
%either a single char, or a cell of char. See AminoAcidProp.csv for the
%tranlation between amino acid letter and its property code.
%
%  PropCode = convAA2PropMEX(AAseq)
%
%  INPUT
%    AAseq: amino acid char
% 
%  OUTPUT
%    PropCode: property code of amino acid char
%
%  EXAMPLE
%    AAseq = 'ARNDCQEGHILKMFPSTWYVZ';
%    PropCode = convAA2PropMEX(AAseq)
%    PropCode = 
%       'HBNASNAGBHHBSFPOOWYHX'
%
%
