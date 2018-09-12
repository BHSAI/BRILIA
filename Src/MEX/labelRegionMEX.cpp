/*  
labelRegionMEX will get a 1xN or Mx1 vector and label form 1 to N the 
number of consecutive-matched numbers.
    
  Label = labelRegionMEX(Region)

  INPUT
    Region: 1xN or Mx1 vector like [0 0 1 1 2 2 2 2 0 0 1 1]
    
  OUTPUT
    Label: vector of size(Region) of regions 1 to N

  EXAMPLE
    Region = [0 0 1 1 2 2 2 2 0 0 1 1];
    Label  = labelRegionMEX(Region)
    Label  = 
             [1 1 2 2 3 3 3 3 4 4 5 5]
*/

#include "mex.h"

void labelDoubleRegion(double *pMat, mwSize Len, double *pLabel) {
    int k = 0;
    double Num = 1;
    for (mwSize j = 0; j < Len; j++) {
        if (pMat[j] == pMat[k]) {
            pLabel[j] = Num;
        } else {
            pLabel[j] = ++Num;
            k = j;
        }
    }
}

void labelBoolRegion(bool *pMat, mwSize Len, double *pLabel) {
    int k = 0;
    double Num = 1;
    for (mwSize j = 0; j < Len; j++) {
        if (pMat[j] == pMat[k]) {
            pLabel[j] = Num;
        } else {
            pLabel[j] = ++Num;
            k = j;
        }
    }
}

void mexFunction(int nlhs,        mxArray *plhs[],
                 int nrhs, const  mxArray *prhs[]) {

    if (nrhs != 1) {
        mexErrMsgIdAndTxt("mxFindConsecOnes:nrhs", "Incorrect number of inputs. Expected 1.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("mxFindConsecOnes:nlhs", "Incorrect number of outputs. Expected 1.");
    }
    if (mxGetN(prhs[0]) > 1 && mxGetM(prhs[0]) > 1) {
        mexErrMsgIdAndTxt("mxFindConsecOnes:prhs", "Input must be a 1xN or Mx1 vector.");
    }

    mwSize M = mxGetM(prhs[0]);
    mwSize N = mxGetN(prhs[0]);
    mwSize Len = M > N ? M : N;
    plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    
    if (M == 0 || N == 0) { return; }

    double *pLabel = mxGetPr(plhs[0]);
    if (mxIsDouble(prhs[0])) {
        double *pMat = mxGetPr(prhs[0]);
        labelDoubleRegion(pMat, Len, pLabel);
    } else if (mxIsLogical(prhs[0])) {
        bool *pMat = mxGetLogicals(prhs[0]);
        labelBoolRegion(pMat, Len, pLabel);
    }
}