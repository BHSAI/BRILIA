/*  
calcPairDistMEX computes the hamming distance and number of valid hotspot 
mutations from SeqA to SeqB, and vice versa. A valid hotspot mutation is
one that agrees with the motif (WRC, GYW, WA, TW) AND has a pairwise 
substitution pattern that agrees with deaminase subsitution patterns. 
A valid Motif is +1, and a valid pairwise mutation is +1, for a maximum 
point of 2 for classical hotspot mutation. Worse is -2 for disagreeing 
motif AND pairwise mutation. If [valid mut] + [invalid mut] < 0, then SeqA 
cannot be parent of SeqB. It is possible SeqA->SeqB, and SeqB->SeqA return
negative ShmDir scores, as in SeqA and SeqB are not linkable under the same
clonotype. Consecutive mismatches add M^2 score, decrease chance SeqA and 
SeqB are related.

    N  A  C  G  T   Interpretted as:
  ---------------   --------------
N | 0  0  0  0  0   N -> - A C G T
A | 0  0  0  1  1   A -> N - C G T
C | 0 -1  0 -1  1   C -> N A - G T 
G | 0  1 -1  0 -1   G -> N A C - T
T | 0 -1  0 -1  0   T -> N A C G -
 
  [Ham, Motif, Mut, Penalty, ShmDist] = calcPairDistMEX(SeqA, SeqB);
 
  [Ham, Motif, Mut, Penalty, ShmDist] = calcPairDistMEX(SeqList);

  INPUT
    SeqA: nt sequence string
    SeqB: nt sequence string
    SeqList: Mx1 cell of nt sequence strings

  OUTPUT
    Ham: MxM hamming distance matrix between SeqA and SeqB
    Motif: MxM total count of hotspot (+1) and non-hotspot motifs (-1)
    Mut: MxM total count of valid (+1) and invalid (-1) nt subsitutions
    Penalty: MxM of total of Mi^2, where Mi = ith mismatch segment
    ShmDist: MxM BRILIA's SHM distance. 

  NOTE
    Each row is the parent, each column is the child.
 
    The "SHM Distance" is calculated as:
      ShmDist = Ham - (Motif + Mut - Penalty)/4;
 
  EXAMPLE

    Seq{1} = 'AACAACGAACGAA';
    Seq{2} = 'AATAATTAATGAA';
    [Ham, Motif, Mut, Penalty, ShmDist] = calcPairDistMEX(Seq);
    Ham =
         0     4
         4     0
    Motif =
         0     2
         2     0
    Mut =
         0     2
        -1     0
    Penalty =
         0     4
         4     0
    ShmDist =
             0    4.0000
        4.7500         0
 */

#include "mex.h"
#include "HotspotTool.hpp"

void mexFunction(int nlhs,        mxArray *plhs[],
                 int nrhs, const  mxArray *prhs[]) {
    
    if (nrhs != 1 || !mxIsCell(prhs[0])) {
        mexErrMsgIdAndTxt("calcPairDistMEX:input", "Input: Must be Mx1 cell of sequences.");
    }
    
    mwSize const *pDim = mxGetDimensions(prhs[0]);
    mwSize NumSeq = pDim[0]*pDim[1];
    
    double *pHamDist, *pValidMotif, *pValidMut, *pPenalty, *pShmDist;
    
    if (nlhs >= 0) { //Hamming Distance matrix
        plhs[0] = mxCreateDoubleMatrix(NumSeq, NumSeq, mxREAL);
        pHamDist = mxGetPr(plhs[0]);
    }

    if (nlhs >= 2) { //Row-to-Col Valid hotspot motif count
        plhs[1] = mxCreateDoubleMatrix(NumSeq, NumSeq, mxREAL);
        pValidMotif = mxGetPr(plhs[1]);
    }
    
    if (nlhs >= 3) { //Row-to-Col Valid pairwise mutation count
        plhs[2] = mxCreateDoubleMatrix(NumSeq, NumSeq, mxREAL);
        pValidMut = mxGetPr(plhs[2]);
    }

    if (nlhs >= 4) { //Penalty from consec mismatch)
        plhs[3] = mxCreateDoubleMatrix(NumSeq, NumSeq, mxREAL);
        pPenalty = mxGetPr(plhs[3]);
    }

    if (nlhs >= 5) { //Compute the SHM distance
        plhs[4] = mxCreateDoubleMatrix(NumSeq, NumSeq, mxREAL);
        pShmDist = mxGetPr(plhs[4]);
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

            mwSize Len = LenA < LenB ? LenA : LenB;
            int Idx1 = j + k*NumSeq;
            int Idx2 = k + j*NumSeq;
            double *pScore = calcSeqShmScore(pSeqA, pSeqB, Len);
            pHamDist[Idx1] = pScore[0];
            pHamDist[Idx2] = pScore[0];
            
            if (nlhs >= 2) {
                pValidMotif[Idx1] = pScore[1];
                pValidMotif[Idx2] = pScore[2];
                
                if (nlhs >= 3) {
                    pValidMut[Idx1] = pScore[3];
                    pValidMut[Idx2] = pScore[4];
                    
                    if (nlhs >= 4) {
                        pPenalty[Idx1] = pScore[5];
                        pPenalty[Idx2] = pScore[5];
                        
                        if (nlhs >= 5) {
                            pShmDist[Idx1] = pScore[0] - (pScore[1] + pScore[3] - pScore[5])/4;
                            pShmDist[Idx2] = pScore[0] - (pScore[2] + pScore[4] - pScore[5])/4;
                        }
                    }
                }
            }
        }
    }
}