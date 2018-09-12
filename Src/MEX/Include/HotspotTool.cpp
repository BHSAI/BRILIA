/*
 *  HotspotTool.cpp contains all the source codes for working with SHM 
 *  hotspots detection for a single sequence strand. 
 */ 

#include "HotspotTool.hpp"
#include "SeqTool.hpp"  // Need to fix compileMex.m to be able to search header file dependencies.
#include "AlignTool.hpp"

//row is NACGT initial letter
//col is NACGT final letter   
double pSHM_TENDENCY[5][5] = {
    { 0,  0,  0,  0,  0}, //N -> - A C G T
    { 0,  0,  0,  1,  1}, //A -> N - C G T
    { 0, -1,  0, -1,  1}, //C -> N A - G T 
    { 0,  1, -1,  0, -1}, //G -> N A C - T
    { 0, -1,  0, -1,  0}  //T -> N A C G -
}; 

bool isWRC(mxChar *pSeq, mwSize Len, mwSize Pos) {
    if (pSeq[Pos] == 'C') {
        if (Pos < 2) { return false; } //Cannot determine
        if (Pos >= 1 && !(pSeq[Pos-1] == 'A' || pSeq[Pos-1] == 'G'))                       { return false; } //check 1st left pos
        if (Pos >= 2 && !(pSeq[Pos-2] == 'A' || pSeq[Pos-2] == 'T' || pSeq[Pos-2] == 'U')) { return false; } //check 2nd left pos
        return true;
    } else {
        return false;
    }
}

bool isGYW(mxChar *pSeq, mwSize Len, mwSize Pos) {
    if (pSeq[Pos] == 'G') {
        if (Pos >= Len-2) { return false; } //Cannot determine
        if (Pos < Len-1 && !(pSeq[Pos+1] == 'C' || pSeq[Pos+1] == 'T' || pSeq[Pos+1] == 'U')) { return false; } //check 1st right pos
        if (Pos < Len-2 && !(pSeq[Pos+2] == 'A' || pSeq[Pos+2] == 'T' || pSeq[Pos+2] == 'U')) { return false; } //check 2nd right pos
        return true;
    } else {
        return false;
    }
}

bool isWA(mxChar *pSeq, mwSize Len, mwSize Pos) {
    if (pSeq[Pos] == 'A') {
        if (Pos < 1) { return false; } //Cannot determine
        if (Pos >= 1 && !(pSeq[Pos-1] == 'A' || pSeq[Pos-1] == 'T' || pSeq[Pos-1] == 'U')) { return false; } //check 1nd left pos
        return true;
    } else {
        return false;
    }
}

bool isTW(mxChar *pSeq, mwSize Len, mwSize Pos) {
    if (pSeq[Pos] == 'T' || pSeq[Pos] == 'U') {
        if (Pos >= Len-1) { return false; } //Cannot determine
        if (Pos < Len-1 && !(pSeq[Pos+1] == 'A' || pSeq[Pos+1] == 'T' || pSeq[Pos+1] == 'U')) { return false; } //check 2nd right pos
        return true;
    } else {
        return false;
    }
}

bool isHotspot(mxChar *pSeq, mwSize Len, mwSize Pos) {
    switch (pSeq[Pos]) {
        case 'A': return isWA(pSeq, Len, Pos);
        case 'C': return isWRC(pSeq, Len, Pos);
        case 'G': return isGYW(pSeq, Len, Pos);
        case 'T': return isTW(pSeq, Len, Pos);
        case 'U': return isTW(pSeq, Len, Pos);
        case 'N': return false;
        case 'X': return false;
        default: mexErrMsgIdAndTxt("HotspotTool_isHotspot:input", "Unknown character detected. Only ACGTUNX are allowed. No lower case.");
    }
}

double labelHotspot(mxChar *pSeq, mwSize Len, mwSize Pos) {
    switch (pSeq[Pos]) {
        case 'A': return isWA(pSeq, Len, Pos)  ? 1 : 0;
        case 'C': return isWRC(pSeq, Len, Pos) ? 2 : 0;
        case 'G': return isGYW(pSeq, Len, Pos) ? 3 : 0;
        case 'T': return isTW(pSeq, Len, Pos)  ? 4 : 0;
        case 'U': return isTW(pSeq, Len, Pos)  ? 4 : 0;
        case 'N': return 0;
        case 'X': return 0;
        default: mexErrMsgIdAndTxt("HotspotTool_labelHotspot:input", "Unknown character detected. Only ACGTUNX are allowed. No lower case.");
    }
}

double countHotspots(mxChar *pSeq, mwSize Len) {
    double Count = 0;
    for (mwSize j = 0; j < Len; j++) {
        Count += isHotspot(pSeq, Len, j) ? 1 : 0;
    }
    return Count;
}

double countHotspots(mxChar *pSeq, mwSize Len, double *pLabel) {
    double Count = 0;
    for (mwSize j = 0; j < Len; j++) {
        pLabel[j] = labelHotspot(pSeq, Len, j);
        Count += pLabel[j] > 0 ? 1 : 0;
    }
    return Count;
}

void calcSeqShmScore(mxChar *pSeqA, mxChar *pSeqB, mwSize Len, double *pScore) {
//     Score[0] = 0; //Hamming dist
//     Score[1] = 0; //A to B hotspot motif count, sum(Motifs), where Motifs = 1 for hotspot, -1 otherwise
//     Score[2] = 0; //B to A hotspot motif count, sum(Motifs), where Motifs = 1 for hotspot, -1 otherwise
//     Score[3] = 0; //A to B pairwise motif count, sum(SHM_TENDENCY(A to B))
//     Score[4] = 0; //B to A pairwise motif count, sum(SHM_TENDENCY(B to A))
//     Score[5] = 0; //Penalty, sum(Mi^2) where Mi = length of ith consec. mismatch segment 
    double ConsecM = 0;
    bool pMatch[Len] = {false};
    cmprSeq(pSeqA, pSeqB, Len, 'n', pMatch);
    for (mwSize j = 0; j < Len; j++) {
        if (!pMatch[j]) {
            int A = nt2int(pSeqA[j]);
            int B = nt2int(pSeqB[j]);
            pScore[0] += 1;
            pScore[1] += isHotspot(pSeqA, Len, j) ? 1 : -1;
            pScore[2] += isHotspot(pSeqB, Len, j) ? 1 : -1;
            pScore[3] += pSHM_TENDENCY[A][B];
            pScore[4] += pSHM_TENDENCY[B][A];
            ConsecM += 1;
        } else {
            if (ConsecM > 1) {//Only if you have > 1 consec mismatch
                pScore[5] += (ConsecM * ConsecM) - 1; // subtract 1 because that's the hamming distance. double penalty. 
            }
            ConsecM = 0;
        }
    }
    if (ConsecM > 1) {
        pScore[5] += (ConsecM * ConsecM) - 1;
    }
}