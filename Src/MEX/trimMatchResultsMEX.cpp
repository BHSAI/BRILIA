/*  
trimMatchResultsMEX will take a binary sequence match results and trim
matches (1's) near left, right, both, or no edges until it encounters a 
a block where 3 out of 4 are matches (ex: 1 1 0 1). Trimmed matches are 
0's.
 
  MatchResults = trimMatchResultsMEX(MatchResults, TrimSide)

  INPUT
    MatchResults: 1xN logical array
    TrimSide ['n', 'l', 'r', 'b']: direction to trim. 'n' does nothing.
 
  OUTPUT
    MatchResults: 1xN logical array with 1's turned to 0's where trimmed

  EXAMPLE
    M = [1 0 0 1 1 1 1 0 0 1] > 0;
    TrimSide = 'b';
    M = trimMatchResultsMEX(M, TrimSide)
    M =
        0   0   0   1   1   1   1   0   0   0
*/

#include "AlignTool.hpp"

void mexFunction(int nlhs,        mxArray *plhs[],
                 int nrhs,  const mxArray *prhs[]) {

    if (nrhs < 1 || nrhs > 2) {
        mexErrMsgIdAndTxt("trimMatchResults:nrhs", "Incorrect number of inputs. Min is 1. Max is 2.");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("trimMatchResults:nlhs", "Too many outputs. Max is 1.");
    }
    if (mxIsComplex(prhs[0]) || mxGetM(prhs[0]) > 1 || !mxIsLogical(prhs[0])) {
        mexErrMsgIdAndTxt("trimMatchResults:prhs", "Input1: Match must be a logical 1xN array.");
    }
    if (!mxIsChar(prhs[1]) || mxGetN(prhs[1])*mxGetM(prhs[1]) > 1) {
        mexErrMsgIdAndTxt("trimMatchResults:prhs", "Input2: TrimSide must be a char 'n', 'l', 'r', or 'b'.");
    }

    mwSize Len = mxGetN(prhs[0]);
    mxChar TrimSide = nrhs > 1 ? *mxGetChars(prhs[1]) : 'n';

    plhs[0] = mxDuplicateArray(prhs[0]);
    bool *pMatch = mxGetLogicals(plhs[0]);

    trimMatchResults(pMatch, Len, TrimSide);
}