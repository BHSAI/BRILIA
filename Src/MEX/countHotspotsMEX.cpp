/*  
countHotspotsMEX will count the number of hotspot motifs in a sequence. 
Partial motifs at the sequence edge will NOT be considered a hotspot.

  N = countHotspotsMEX(Seq)
 
  INPUT
    Seq: nt sequence string or cell of sequence strings (all CAPS)

  OUTPUT
    N: scalar or matrix of number of hotspots

  EXAMPLE

    Seq = {'TGCTGCTGC';
           'TAGGUAG'};
    [N, S] = countHotspotsMEX(Seq);

*/

#include "HotspotTool.hpp"

void mexFunction(int nlhs,        mxArray *plhs[],
                 int nrhs, const  mxArray *prhs[]) {
    
    if (nrhs != 1) {
        mexErrMsgIdAndTxt("countHotspotsMEX:nrhs", "Incorrect number of inputs - expected 1 input.");
    }
    if (nlhs < 1 || nlhs > 2) {
        mexErrMsgIdAndTxt("countHotspotsMEX:nlhs", "Incorrect number of outputs. Expected 1 or 2.");
    }
    if (nrhs == 1 && !(mxIsCell(prhs[0]) || mxIsChar(prhs[0]))) {
        mexErrMsgIdAndTxt("countHotspotsMEX:prhs", "Input must be a string or cell of sequence strings.");
    } 
    
    if (mxIsCell(prhs[0])) {
        mwSize const *pDim = mxGetDimensions(prhs[0]);
        mwSize NumSeq = pDim[0]*pDim[1];
        
        plhs[0] = mxCreateDoubleMatrix(pDim[0], pDim[1], mxREAL);
        double *pCount = mxGetPr(plhs[0]);
        if (nlhs == 2) {
            plhs[1] = mxCreateCellMatrix(pDim[0], pDim[1]);
        }
        
        for (int j = 0; j < NumSeq; j++) {
            mxArray *pCell = mxGetCell(prhs[0], j);
            if (!mxIsChar(pCell)) { 
                mexErrMsgIdAndTxt("countHotspotsMEX:input", "Input: Not all cells are strings.");
            }
            mxChar *pSeq = mxGetChars(pCell);
            mwSize Len = mxGetN(pCell);
            if (nlhs == 1) {
                pCount[j] = countHotspots(pSeq, Len);
            } else {
                mxSetCell(plhs[1], j, mxCreateDoubleMatrix(1, Len, mxREAL));
                pCount[j] = countHotspots(pSeq, Len, mxGetPr(mxGetCell(plhs[1], j)));
            }
        }
    } else {
        mxChar *pSeq = mxGetChars(prhs[0]);
        mwSize Len = mxGetN(prhs[0]);
        if (nlhs == 1) {
            plhs[0] = mxCreateDoubleScalar(countHotspots(pSeq, Len));
        } else {
            plhs[1] = mxCreateDoubleMatrix(1, Len, mxREAL);
            plhs[0] = mxCreateDoubleScalar(countHotspots(pSeq, Len, mxGetPr(plhs[1])));
        }
    }
}