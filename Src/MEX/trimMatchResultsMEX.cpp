/*  trimMatchResultsMEX will take a binary sequence match results and trim
 *  matches (1's) near left, right, both, or no edges until it encounters a 
 *  a block where 3 out of 4 are matches (ex: 1 1 0 1). Trimmed matches are 
 *  0's.
 *
 *  MatchResults = trimMatchResultsMEX(MatchResults, TrimSide)
 *
 *  INPUT
 *    MatchResults: 1xN logical array
 *    TrimSide ['l', 'r', 'b', 'n']: direction to trim. 'n' does nothing.
 *
 *  OUTPUT
 *    MatchResults: 1xN logical array with 1's turned to 0's where trimmed
 *
 *  EXAMPLE
 *    M = [1 0 0 1 1 1 1 0 0 1] > 0;
 *    TrimSide = 'b';
 *    M = trimMatchResultsMEX(M, TrimSide)
 *      =  0   0   0   1   1   1   1   0   0   0
 */

#include "mex.h"
#include <string>

void trimMatchResultsMEX(bool *pMatchI,mwSize N,  mxChar TrimSide, bool *pMatchO) {
    if (N <= 4) return;
    double Sum;
    mwSize L(0), R(N-1);

    //Determine left start point
    if (toupper(TrimSide) == 'L' || toupper(TrimSide) == 'B') {
        Sum = 0;
        mwSize i;
        for (i = 0; i < 4; i++) {//Initial 1st 4 sum
            if (pMatchI[i]) Sum++;
        }
        while (i < N) {
            if (Sum >= 3) {
                L = i-3;
                break;
            }
            i++;
            if (pMatchI[i-4]) Sum--;
            if (pMatchI[i]) Sum++;
        }
        for (i = 0; i < L; i++) {pMatchO[i] = 0;}
    }

    //Determine right end point
    if (toupper(TrimSide) == 'R' || toupper(TrimSide) == 'B') {
        Sum = 0;
        mwSize i; //NOTE: mwSize is unsigned, so i will NEVER be negative.
        for (i = N; i >= N-3; i--) {//Initial 1st 4 sum 
            if (pMatchI[i-1]) Sum++;
        }
        while (i >= 1) {
            if (Sum >= 3) {
                R = i+4;
                break;
            }
            i--;
            if (pMatchI[i+3]) Sum--;
            if (pMatchI[i-1]) Sum++;
        }
        for (i = N; i >= R; i--) {pMatchO[i-1] = 0;}
    }
}

void mexFunction(int nlhs,        mxArray *plhs[],
                 int nrhs,  const mxArray *prhs[]) {

    if (nrhs != 2) {
        mexErrMsgIdAndTxt("trimMatchResults:nrhs", "Need 2 inputs.");
    }
    
    if (mxIsComplex(prhs[0]) || mxGetM(prhs[0]) > 1 || !mxIsLogical(prhs[0])) {
        mexErrMsgIdAndTxt("trimMatchResults:prhs", "1st input [Match] must be a logical 1xN vector.");
    }
    
    if (!mxIsChar(prhs[1]) || mxGetN(prhs[1])*mxGetM(prhs[1]) > 1) {
        mexErrMsgIdAndTxt("trimMatchResults:prhs", "2nd input [TrimSide] must be a char 'l', 'r', or 'b'.");
    }

    if (nlhs != 1) {
        mexErrMsgIdAndTxt("trimMatchResults:nlhs", "Need 1 output.");
    }

    mwSize N = mxGetN(prhs[0]);
    bool *pMatchI = mxGetLogicals(prhs[0]);
    mxChar TrimSide = *mxGetChars(prhs[1]);

    plhs[0] = mxDuplicateArray(prhs[0]);
    bool *pMatchO = mxGetLogicals(plhs[0]);

    trimMatchResultsMEX(pMatchI, N, TrimSide, pMatchO);
}

