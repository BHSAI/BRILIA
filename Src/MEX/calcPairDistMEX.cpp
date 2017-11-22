/*  
calcShmHamDistMEX computes the hamming distance and SHM distance from seqA 
to seqB, and vice versa. 

Score = calcShmHamDistMEX(SeqA, SeqB)

INPUT

OUTPUT

EXAMPLE
*/

#include "mex.h"
#include <string>

int nt2int(mxChar charA) {
    switch (charA) {
        case 'n':
            return 0;
        case 'N':
            return 0;
        case 'x':
            return 0;
        case 'X':
            return 0;
        case 'a':
            return 1;
        case 'A':
            return 1;
        case 'c':
            return 2;
        case 'C':
            return 2;
        case 'g':
            return 3;
        case 'G':
            return 3;
        case 't':
            return 4;
        case 'T':
            return 4;
        default:
            return 5;
    }
}

void calcShmHamDistMEX(mxChar *pSeqA, mxChar *pSeqB, mwSize MinLen, int pShmTendency[5][5], int nlhs, double *pHamDist, double *pShmDistAtoB, double *pShmDistBtoA) {
    double Penalty = 0;
    double Misses  = 0;
    double AtoBTendency = 0;
    double BtoATendency = 0;
    int numA, numB;
    
    if (nlhs > 1) {
        for (int j = 0; j < MinLen; j++) {
            numA = nt2int(pSeqA[j]);
            numB = nt2int(pSeqB[j]);
            if (numA == 0 || numB == 0) { continue; }
            if (numA == numB) {
                Penalty += Misses*Misses;
                Misses = 0;
            } else {
                *pHamDist = *pHamDist+1;
                Misses++;
            }
            AtoBTendency += pShmTendency[numA][numB];
            BtoATendency += pShmTendency[numB][numA];
        }
        Penalty += Misses*Misses;
        
        *pShmDistAtoB = Penalty - 0.5*AtoBTendency;
        *pShmDistBtoA = Penalty - 0.5*BtoATendency;
    } else {
        for (int j = 0; j < MinLen; j++) {
            numA = nt2int(pSeqA[j]);
            numB = nt2int(pSeqB[j]);
            if (numA ==0 || numB == 0) { continue; }
            if (numA != numB) {
                *pHamDist = *pHamDist+1;
            }
        }
    }
}

void mexFunction(int nlhs,        mxArray *plhs[],
                 int nrhs, const  mxArray *prhs[]) {
    
    if (nrhs < 1 || nrhs > 2) {
        mexErrMsgIdAndTxt("calcShmHamDistMEX:input", "Need 1 to 2 inputs.");
    }
    
    bool IsCell = false;
    if (nrhs == 1) {
        if (!mxIsCell(prhs[0])) {
            mexErrMsgIdAndTxt("calcShmHamDistMEX:input", "Input1: Must be a cell array of strings.");
        } else {
            IsCell = true;
            if (nlhs > 2)
                mexErrMsgIdAndTxt("calcShmHamDistMEX:output", "No 3rd output when using cell array inputs of sequences.");
        }
    } else {
        IsCell = false;
        if (!mxIsChar(prhs[0]))
            mexErrMsgIdAndTxt("calcShmHamDistMEX:input", "Input1: SeqA must be a string.");
        if (!mxIsChar(prhs[1]))
            mexErrMsgIdAndTxt("calcShmHamDistMEX:input", "Input2: SeqB must be a string.");
        if (nlhs > 3)
            mexErrMsgIdAndTxt("calcShmHamDistMEX:input", "No 4th output when using string array inputs of sequences.");
    }

    mwSize const *pDim = mxGetDimensions(prhs[0]);
    mwSize NumSeq = pDim[0]*pDim[1];
                
    //row is NACGT initial letter
    //col is NACGT final letter
    int ShmTendency[5][5] = {
        { 0,  0,  0,  0,  0},
        { 0,  0,  0,  1,  1},
        { 0, -1,  0, -1,  1},
        { 0,  1, -1,  0, -1},
        { 0, -1,  0, -1,  0}
    }; 
    
    if (IsCell) {
        double *pHamDist, *pShmDist;
        if (nlhs >= 1) { //Hamming Distance matrix
            plhs[0] = mxCreateDoubleMatrix(NumSeq, NumSeq, mxREAL);
            pHamDist = mxGetPr(plhs[0]);
        }

        if (nlhs >= 2) { //A to B distance on lower left. B to A distance on upper right.
            plhs[1] = mxCreateDoubleMatrix(NumSeq, NumSeq, mxREAL);
            pShmDist = mxGetPr(plhs[1]);
        }

        for (int j = 0; j < NumSeq; j++) {
            mxArray *pCellA = mxGetCell(prhs[0], j);
            if (!mxIsChar(pCellA)) { continue; }
            mxChar *pSeqA = mxGetChars(pCellA);
            mwSize LenA = mxGetN(pCellA);
            
            for (int k = 0; k < j; k++) {
                mxArray *pCellB = mxGetCell(prhs[0], k);
                if (!mxIsChar(pCellB)) { continue; }
                mxChar *pSeqB = mxGetChars(pCellB);
                mwSize LenB = mxGetN(pCellB);

                mwSize MinLen = LenA < LenB ? LenA : LenB;
                int Idx1 = j + k*NumSeq;
                int Idx2 = k + j*NumSeq;
                calcShmHamDistMEX(pSeqA, pSeqB, MinLen, ShmTendency, nlhs, &pHamDist[Idx1], &pShmDist[Idx1], &pShmDist[Idx2]);
                pHamDist[Idx2] = pHamDist[Idx1];
            }
        }
    } else {
        double HamDist = 0, AtoBDist = 0, BtoADist = 0;        
        mxChar *pSeqA = mxGetChars(prhs[0]);
        mwSize LenA = mxGetN(prhs[0]);
        mxChar *pSeqB = mxGetChars(prhs[1]);
        mwSize LenB = mxGetN(prhs[1]);
        mwSize MinLen = LenA < LenB ? LenA : LenB;
        calcShmHamDistMEX(pSeqA, pSeqB, MinLen, ShmTendency, nlhs, &HamDist, &AtoBDist, &BtoADist);

        if (nlhs >= 1)  //Hamming Distance
            plhs[0] = mxCreateDoubleScalar(HamDist);
        if (nlhs >= 2)  //SeqA as parent, SeqB as child (Par2Child)
            plhs[1] = mxCreateDoubleScalar(AtoBDist);
        if (nlhs >= 3)  //SeqB as parent, SeqA as child (Child2Par)
            plhs[2] = mxCreateDoubleScalar(BtoADist);
    }
}