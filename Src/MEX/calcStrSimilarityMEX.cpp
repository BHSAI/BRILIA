/*
 calcStrSimilarityMEX will compare two strings to determine the % 
 similarity with respect to the shorter sequences.
 
 Example:
   A = 'ACGTTACGTT';
   B = 'ACGT';
   [c, d] = calcStrSimilarityMEX(A, B)
    
   C = {'ACGT', 'ACTTGCA', 'ACGTTACGTTCA'}
   [c, d] = calcStrSimilarityMEX(C)
  
   
 */

#include "mex.h"
#include <string>

void calcStrHammingMEX(mxChar *pSeqA, mxChar *pSeqB, mwSize LenA, mwSize LenB, double *pHamDist) {
    mwSize MinLen = LenA < LenB ? LenA : LenB;
    for (int j = 0; j < MinLen; j++) {
        if (tolower(pSeqB[j]) == tolower(pSeqA[j]))
            *pHamDist = *pHamDist + 1;
    }
}

void calcStrSimilarityMEX(mxChar *pSeqA, mxChar *pSeqB, mwSize LenA, mwSize LenB, double *pHamDist, double *pSimDist) {
    mwSize MinLen = LenA < LenB ? LenA : LenB;
    for (int j = 0; j < MinLen; j++) {
        if (tolower(pSeqB[j]) == tolower(pSeqA[j]))
            *pHamDist = *pHamDist + 1;
    }
    *pSimDist = *pHamDist / (double) MinLen;
}

void mexFunction(int nlhs,        mxArray *plhs[],
                 int nrhs, const  mxArray *prhs[]) {
    
    if (nrhs < 1 || nrhs > 2) {
        mexErrMsgIdAndTxt("calcStrSimilarityMEX:input", "Need 1 to 2 inputs.");
    }
    
    if (nlhs > 2) {
        mexErrMsgIdAndTxt("calcStrSimilarityMEX:input", "Too many outputs (max 2).");
    }
    
    bool IsCell = false;
    if (nrhs == 1) {
        if (!mxIsCell(prhs[0])) {
            mexErrMsgIdAndTxt("calcStrSimilarityMEX:input", "Input1: Must be a cell array of strings.");
        } else {
            IsCell = true;
        }
    } else {
        IsCell = false;
        if (!mxIsChar(prhs[0]))
            mexErrMsgIdAndTxt("calcStrSimilarityMEX:input", "Input1: SeqA must be a string.");
        if (!mxIsChar(prhs[1]))
            mexErrMsgIdAndTxt("calcStrSimilarityMEX:input", "Input2: SeqB must be a string.");
    }

    if (IsCell) {
        mwSize const *pDim = mxGetDimensions(prhs[0]);
        mwSize NumSeq = pDim[0]*pDim[1];
        
        if (nlhs == 1) {
            plhs[0] = mxCreateDoubleMatrix(NumSeq, NumSeq, mxREAL);
            double *pHamDist = mxGetPr(plhs[0]);
            for (int j = 0; j < NumSeq; j++) {
                mxArray *pCellA = mxGetCell(prhs[0], j);
                if (!mxIsChar(pCellA)) { continue; }
                mxChar *pSeqA = mxGetChars(pCellA);
                mwSize   LenA = mxGetN(pCellA);
                for (int k = 0; k < j; k++) {
                    mxArray *pCellB = mxGetCell(prhs[0], k);
                    if (!mxIsChar(pCellB)) { continue; }
                    mxChar *pSeqB = mxGetChars(pCellB);
                    mwSize   LenB = mxGetN(pCellB);
                    int Idx1 = j + k*NumSeq;
                    int Idx2 = k + j*NumSeq;
                    calcStrHammingMEX(pSeqA, pSeqB, LenA, LenB, &pHamDist[Idx1]);
                    pHamDist[Idx2] = pHamDist[Idx1];
                }
            }
        }
        
        if (nlhs == 2) {
            plhs[0] = mxCreateDoubleMatrix(NumSeq, NumSeq, mxREAL);
            double *pHamDist = mxGetPr(plhs[0]);
            plhs[1] = mxCreateDoubleMatrix(NumSeq, NumSeq, mxREAL);
            double *pSimDist = mxGetPr(plhs[1]);

            for (int j = 0; j < NumSeq; j++) {
                mxArray *pCellA = mxGetCell(prhs[0], j);
                if (!mxIsChar(pCellA)) { continue; }
                mxChar *pSeqA = mxGetChars(pCellA);
                mwSize   LenA = mxGetN(pCellA);
                for (int k = 0; k < j; k++) {
                    mxArray *pCellB = mxGetCell(prhs[0], k);
                    if (!mxIsChar(pCellB)) { continue; }
                    mxChar *pSeqB = mxGetChars(pCellB);
                    mwSize   LenB = mxGetN(pCellB);
                    int Idx1 = j + k*NumSeq;
                    int Idx2 = k + j*NumSeq;
                    calcStrSimilarityMEX(pSeqA, pSeqB, LenA, LenB, &pHamDist[Idx1], &pSimDist[Idx1]);
                    pHamDist[Idx2] = pHamDist[Idx1];
                    pSimDist[Idx2] = pSimDist[Idx1];
                }
            }
        }
        
    } else {
        double HamDist = 0, SimDist = 0;
        mxChar *pSeqA = mxGetChars(prhs[0]);
        mwSize   LenA = mxGetN(prhs[0]);
        mxChar *pSeqB = mxGetChars(prhs[1]);
        mwSize   LenB = mxGetN(prhs[1]);
        calcStrSimilarityMEX(pSeqA, pSeqB, LenA, LenB, &HamDist, &SimDist);

        if (nlhs >= 1)
            plhs[0] = mxCreateDoubleScalar(HamDist);
        if (nlhs >= 2)
            plhs[1] = mxCreateDoubleScalar(SimDist);
    }
}