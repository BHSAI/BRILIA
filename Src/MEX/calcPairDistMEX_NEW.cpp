/*  
calcPairDistMEX computes the hamming distance and SHM distance from SeqA 
to SeqB, and vice versa. 

  [Ham, A2B, B2A] = calcShmHamDistMEX(SeqA, SeqB)

  [HamDist, ShmDist] = calcShmHamDistMEX(SeqList)
 
  INPUT
    SeqA: nt sequence string
    SeqB: nt sequence string
    SeqList: Mx1 cell of nt sequence strings

  OUTPUT
    Ham: hamming distance between SeqA and SeqB
    A2B: BRILIA custom SHM distance from SeqA(parent) to SeqB(child)
    B2A: BRILIA custom SHM distance from SeqB(child) to SeqA(parent)
    HamDist: Matrix of hamming distance between sequences in SeqList
    ShmDist: Matrix of BRILIA's SHM distance. Each row is the parent, 
      each column is the child.

  EXAMPLE
    SeqA = 'ACGGAGATGAACAGT'
    SeqB = 'ATGGAGACGAGTAGT'
    [Ham, P2C, C2P] = calcPairDistMEX(SeqA, SeqB);
    Ham =
          4
    P2C =
          4.5000
    C2P =
          5
 *
 *  [Ham, P2C, C2P] = calcPairDistMEX('AACG', 'AATG')
 
    SeqList = {'ACGGAGATGAACAGT', 'ATGGAGACGAGTAGT', 'ATGGATGCGAGTGGA'}
 
    SeqList = {'AAC' 'AAT'}
    [HamDist, ShmDist] = calcPairDistMEX(SeqList);
    HamDist =
             0      4      8
             4      0      4
             8      4      0

    ShmDist =
           0.0    4.5   18.5
           5.0    0.0    6.0
          18.0    5.0    0.0

*/

#include "mex.h"
#include <string>

//row is NACGT initial letter
//col is NACGT final letter   
double pShmTendency[5][5] = {
    { 0,  0,  0,  0,  0}, //N -> - A C G T
    { 0,  0,  0,  1,  1}, //A -> N - C G T
    { 0, -1,  0, -1,  1}, //C -> N A - G T 
    { 0,  1, -1,  0, -1}, //G -> N A C - T
    { 0, -1,  0, -1,  0}  //T -> N A C G -
}; 

//WARNING: we don't do lower case, ever.
double nt2int(mxChar charA) {
    switch (charA) {
        case 'N': return 0;
        case 'X': return 0;
        case 'A': return 1;
        case 'C': return 2;
        case 'G': return 3;
        case 'T': return 4;
        case 'U': return 4;
        default: mexErrMsgIdAndTxt("calcPairDistMEX:input", "Unknown character detected. Only NXACGTU are allowed. No lower case.");
    }
}

bool isWRC(mxChar *pSeq, mwSize Len, mwSize Pos) {
    if (pSeq[Pos] == 'C') {
        if (Pos >= 1 && !(pSeq[Pos-1] == 'A' || pSeq[Pos-1] == 'G'))                       { return false; } //check 1st left pos
        if (Pos >= 2 && !(pSeq[Pos-2] == 'A' || pSeq[Pos-2] == 'T' || pSeq[Pos-2] == 'U')) { return false; } //check 2nd left pos
        return true;
    } else {
        return false;
    }
}

bool isGYW(mxChar *pSeq, mwSize Len, mwSize Pos) {
    if (pSeq[Pos] == 'G') {
        if (Pos < Len-1 && !(pSeq[Pos+1] == 'C' || pSeq[Pos+1] == 'T' || pSeq[Pos+1] == 'U')) { return false; } //check 1st right pos
        if (Pos < Len-2 && !(pSeq[Pos+2] == 'A' || pSeq[Pos+2] == 'T' || pSeq[Pos+2] == 'U')) { return false; } //check 2nd right pos
        return true;
    } else {
        return false;
    }
}

bool isWA(mxChar *pSeq, mwSize Len, mwSize Pos) {
    if (pSeq[Pos] == 'A') {
        if (Pos >= 1 && !(pSeq[Pos-1] == 'A' || pSeq[Pos-1] == 'T' || pSeq[Pos-1] == 'U')) { return false; } //check 1nd left pos
        return true;
    } else {
        return false;
    }
}

bool isTW(mxChar *pSeq, mwSize Len, mwSize Pos) {
    if (pSeq[Pos] == 'T') {
        if (Pos < Len-1 && !(pSeq[Pos+1] == 'A' || pSeq[Pos+1] == 'T' || pSeq[Pos+1] == 'U')) { return false; } //check 2nd right pos
        return true;
    } else {
        return false;
    }
}

double getHotspotScore(mxChar *pSeq, mwSize Len, mwSize Pos) {
    switch (pSeq[Pos]) {
        case 'C': 
            return isWRC(pSeq, Len, Pos) ? 1 : -1;
        case 'G':
            return isGYW(pSeq, Len, Pos) ? 1 : -1;
        case 'T':
            return isTW(pSeq, Len, Pos)  ? 1 : -1;
        case 'A': 
            return isWA(pSeq, Len, Pos)  ? 1 : -1;
    }
}

// MutType will only determine if if the mutation agrees/disagrees with SHM
//   1 for known hotspot, -1 for unknown hotspot
//   1 for known pairwise mutations, 0 for neutral, and -1 for unfavorable mutations
// Thus, the total sum of scores can be computed as:
//    2 = both hotspot motif AND pairwise mutation are correct
//    1 = correct hotspot BUT neutral pairwise mutation
//    0 = correct hotspot AND incorrect mutation (vice versa), OR cannot determine
//   -1 = incorrect hotspot BUT neutral pairwise mutation
//   -2 = incorrect hotspot AND incorrect mutations
void getMutType(mxChar *pSeqA, mxChar *pSeqB, mwSize Len, mwSize Pos, double *pAtoB, double *pBtoA) {
    int numA = nt2int(pSeqA[Pos]);
    int numB = nt2int(pSeqB[Pos]);
    *pAtoB = pShmTendency[numA][numB] + getHotspotScore(pSeqA, Len, Pos);
    *pBtoA = pShmTendency[numB][numA] + getHotspotScore(pSeqB, Len, Pos);
}

double calcShmHamDist(mxChar *pSeqA, mxChar *pSeqB, mwSize Len, double *pHamDist, double *pShmDistAtoB, double *pShmDistBtoA) {
    double Value = 0;
    double AtoB = 0;
    double BtoA = 0;
    double Penalty = 0;
    for (mwSize j = 0; j < Len; j++) {
        if (pSeqA[j] != pSeqB[j]) {
            getMutType(pSeqA, pSeqB, Len, j, &AtoB, &BtoA);
            *pHamDist += 1;
            *pShmDistAtoB -= AtoB + Penalty;
            *pShmDistBtoA -= BtoA + Penalty;
            Penalty += 4; // increment 4 since you'll divide by 4 later. Consec mismatches add a lot of penalty.
        } else {
            Penalty = 0;
        }
    }
    
    //Final calculation is HAM + 0.25*XtoYdist
    *pShmDistAtoB = *pHamDist + 0.25* *pShmDistAtoB; 
    *pShmDistBtoA = *pHamDist + 0.25* *pShmDistBtoA;
}
/*   
void calcShmHamDistMEX(mxChar *pSeqA, mxChar *pSeqB, mwSize MinLen, int nlhs, double *pHamDist, double *pShmDistAtoB, double *pShmDistBtoA) {
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
            if (numA == 0 || numB == 0) { continue; }
            if (numA != numB) {
                *pHamDist = *pHamDist+1;
            }
        }
    }
}
*/
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
                calcShmHamDist(pSeqA, pSeqB, MinLen, &pHamDist[Idx1], &pShmDist[Idx1], &pShmDist[Idx2]);
                pHamDist[Idx2] = pHamDist[Idx1];
            }
        }
    } else {
        double HamDist = 0, AtoBDist = 0, BtoADist = 0;        
        mxChar *pSeqA = mxGetChars(prhs[0]);
        mxChar *pSeqB = mxGetChars(prhs[1]);
        mwSize LenA = mxGetN(prhs[0]);
        mwSize LenB = mxGetN(prhs[1]);
        mwSize MinLen = LenA < LenB ? LenA : LenB;
        calcShmHamDist(pSeqA, pSeqB, MinLen, &HamDist, &AtoBDist, &BtoADist);

        if (nlhs >= 1)  //Hamming Distance
            plhs[0] = mxCreateDoubleScalar(HamDist);
        if (nlhs >= 2)  //SeqA as parent, SeqB as child (Par2Child)
            plhs[1] = mxCreateDoubleScalar(AtoBDist);
        if (nlhs >= 3)  //SeqB as parent, SeqA as child (Child2Par)
            plhs[2] = mxCreateDoubleScalar(BtoADist);
    }
}