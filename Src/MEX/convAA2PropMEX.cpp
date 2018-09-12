/*
convAA2PropMEX will translate AA seq to its property code. Input is
either a single char, or a cell of char. See AminoAcidProp.csv for the
tranlation between amino acid letter and its property code.

  PropCode = convAA2PropMEX(AAseq)

  INPUT
    AAseq: amino acid char
 
  OUTPUT
    PropCode: property code of amino acid char

  EXAMPLE
    AAseq = 'ARNDCQEGHILKMFPSTWYVZ';
    PropCode = convAA2PropMEX(AAseq)
    PropCode = 
       'HBNASNAGBHHBSFPOOWYHX'
*/

#include "StringTool.hpp"

void convAA2Prop(mxChar *pSeq, mwSize Len, mxChar *pProp) {
    for (mwSize j = 0; j < Len; j++) {
        if        (containsChar(pSeq[j], "AILV")){
            pProp[j] = 'H';
        } else if (containsChar(pSeq[j], "RHK")){
            pProp[j] = 'B';
        } else if (containsChar(pSeq[j], "DE")){
            pProp[j] = 'A';
        } else if (containsChar(pSeq[j], "NQ")){
            pProp[j] = 'N';
        } else if (containsChar(pSeq[j], "ST")){
            pProp[j] = 'O';
        } else if (containsChar(pSeq[j], "CM")){
            pProp[j] = 'S';
        } else if (containsChar(pSeq[j], "PWYGF")){
            pProp[j] = pSeq[j];
        } else {
            pProp[j] = 'X';
        }
    }
}

void mexFunction(int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    
    if (nrhs != 1) {
        mexErrMsgIdAndTxt("convAA2PropMEX:nrhs", "Incorrect number of inputs. Expect 1.");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("convAA2PropMEX:nlhs", "Too many outputs. Max is 1.");
    }
    if (!mxIsChar(prhs[0]) || mxGetM(prhs[0]) > 1) {
        mexErrMsgIdAndTxt("convAA2PropMEX:prhs", "Input must be 1xN char array.");
    }
    
    mxChar *pSeq = mxGetChars(prhs[0]);
    mwSize Len = mxGetN(prhs[0]);
    
    mwSize Dim[2] = {1, Len};
    plhs[0] = mxCreateCharArray(2, Dim);
    mxChar *pProp = mxGetChars(plhs[0]);
    
    convAA2Prop(pSeq, Len, pProp);
}