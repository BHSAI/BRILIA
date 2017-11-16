#include "mex.h"

void mxFindConsecOnes(double *MatI, mwSize M, double *MatO) {
    mwSize i;
    double j = 0;
    bool b = false;
    for (i = 0; i < M; i++) {
        if (!b && MatI[i] > 0) {
            b = true;
            j++;
        } else if (b && MatI[i] == 0) {
            b = false;
        }
        if (MatI[i] == 1) {
            MatO[i] = j;
        }
    }
}

void mxFindConsecOnes(bool *MatI, mwSize M, double *MatO) {
    mwSize i;
    double j = 0;
    bool b = false;
    for (i = 0; i < M; i++) {
        if (!b && MatI[i]) {
            b = true;
            j++;
        } else if (b && !MatI[i]) {
            b = false;
        }
        if (MatI[i]) {
            MatO[i] = j;
        }
    }
}

void mexFunction(int nlhs,        mxArray *plhs[],
                 int nrhs, const  mxArray *prhs[]) {

    if (nrhs != 1) {
        mexErrMsgIdAndTxt("mxFindConsecOnes:nrhs", "Need one input only.");
    }
    
    if (mxIsComplex(prhs[0])) { 
        mexErrMsgIdAndTxt("mxFindConsecOnes:prhs", "Must be a real number.");
    }
    
    if (mxGetN(prhs[0]) > 1 && mxGetM(prhs[0]) > 1) {
        mexErrMsgIdAndTxt("mxFindConsecONes:prhs", "Must be a 1xN or Mx1 vector.");
    }
    
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("mxFindConsecOnes:nlhs", "Need one output only.");
    }

    mwSize M = mxGetM(prhs[0]);
    mwSize N = mxGetN(prhs[0]);
    mwSize P = M > N ? M : N;
    double *MatO;
    plhs[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    MatO = mxGetPr(plhs[0]);

    if (mxIsDouble(prhs[0])) {
        double *MatI; 
        MatI = mxGetPr(prhs[0]);
        mxFindConsecOnes(MatI, P, MatO);
    } else if (mxIsLogical(prhs[0])) {
        bool *MatI;
        MatI = mxGetLogicals(prhs[0]);
        mxFindConsecOnes(MatI, P, MatO);
    }
}

