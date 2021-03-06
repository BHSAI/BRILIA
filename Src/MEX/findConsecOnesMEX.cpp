/*  
findConsecOnesMEX will get a 1xN or Mx1 vector and label form 1 to N the
number of consecutive-matched ones.

  Label = findConsecOnesMEX(Region)
    
  INPUT
    Region: 1xN or Mx1 vector like [0 0 1 1 1 1 0 0 0 0 1 1]
    
  OUTPUT
    Label: vector of size(Region) of regions 1 to N

  EXAMPLE
    Region = [0 0 1 1 1 1 0 0 0 0 1 1];
    Label  = findConsecOnesMEX(Region)
    Label  = 
             [0 0 1 1 1 1 0 0 0 0 2 2]
*/

#include "mex.h"

void findConsecOnes(double *MatI, mwSize M, double *MatO) {
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

void findConsecOnes(bool *MatI, mwSize M, double *MatO) {
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
        mexErrMsgIdAndTxt("mxFindConsecOnes:nrhs", "Incorrect number of inputs. Expected 1.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("mxFindConsecOnes:nlhs", "Incorrect number of inputs. Expected 1.");
    }
    if (mxGetN(prhs[0]) > 1 && mxGetM(prhs[0]) > 1) {
        mexErrMsgIdAndTxt("mxFindConsecONes:prhs", "Input must be a 1xN or Mx1 vector.");
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
        findConsecOnes(MatI, P, MatO);
    } else if (mxIsLogical(prhs[0])) {
        bool *MatI;
        MatI = mxGetLogicals(prhs[0]);
        findConsecOnes(MatI, P, MatO);
    }
}

