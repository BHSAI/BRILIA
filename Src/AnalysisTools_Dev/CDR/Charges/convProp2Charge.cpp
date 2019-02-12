/*  
convProp2Charge takes 2 amino acid property strings and then returns the 
net change in charge. Mainly looking for B(+1), A(-1), and all others (0).

[ChargeArray, TotCharge] = convProp2Charge(Prop1, Prop2)

 INPUT
   Prop1: starting AA property code
   Prop2: ending AA property code
 
 OUTPUT
   ChargeArray: 1xN double matrix of charge differences
   TotCharge: (charges in Prop 2) - (charges in Prop 1)

 EXAMPLE
   Prop1 = 'HAHHAAH' %(-3)
   Prop2 = 'HBHHBBH' %(+3)
   [C, D] = convProp2Charge(Prop1, Prop2)
   C =
       0   2   0   0   2   2   0
   D =
       6
 */

#include "mex.h"
#include <string>

double getCharge(mxChar AA) {
    double Charge = 0;
    switch(AA) {
        case 'A': Charge = -1;
                  break;
        case 'B': Charge =  1;
                  break;
    }
    return Charge;
}

double sumCharges(double *pCharges, mwSize Len) {
    double Score = 0;
    for (mwSize i = 0; i < Len; i++) {
        Score += pCharges[i];
    }
    return Score;
}

void cmprCharges(mxChar *pSeqA, mxChar *pSeqB, mwSize Len, double *pCharges) {
    for (mwSize i = 0; i < Len; i++) {
        pCharges[i] = getCharge(pSeqB[i]) - getCharge(pSeqA[i]);
    }
} 

void mexFunction(int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {

    if (nrhs < 2)
        mexErrMsgIdAndTxt("convProp2ChargeMEX:input", "Not enough inputs.");

    if (nlhs > 2)
        mexErrMsgIdAndTxt("convProp2ChargeMEX:output", "Too many outputs.");
    
    if (!mxIsChar(prhs[0]) || mxGetM(prhs[0]) > 1 || !mxIsChar(prhs[1]) || mxGetM(prhs[1]) > 1)
        mexErrMsgIdAndTxt("convProp2ChargeMEX:input", "Input: AA prop sequence should be Char.");

    if (mxGetN(prhs[0]) != mxGetN(prhs[1])) 
        mexErrMsgIdAndTxt("convProp2ChargeMEX:input", "Input: Both AA sequences should be of same lengths.");

    mxChar *pSeqA = mxGetChars(prhs[0]);
    mxChar *pSeqB = mxGetChars(prhs[1]);
    mwSize Len    = mxGetN(prhs[0]);
    
    plhs[0] = mxCreateDoubleMatrix(1, Len, mxREAL);
    double *pCharges = mxGetPr(plhs[0]);
    cmprCharges(pSeqA, pSeqB, Len, pCharges);
    
    if (nlhs == 2) {
        double Score = sumCharges(pCharges, Len);
        plhs[1] = mxCreateDoubleScalar(Score);
    }
}