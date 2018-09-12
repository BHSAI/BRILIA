/*  SeqTool contains all codes related to NT<->AA, AA<->Int, NT<->AA 
 *  conversions.
 *  WARNING: No lowercase sequence letter is ever used in this tool set!
 */

#include "SeqTool.hpp"

int nt2int(mxChar nt) {
    switch (nt) {
        case 'N': return 0;
        case 'X': return 0;
        case 'A': return 1;
        case 'C': return 2;
        case 'G': return 3;
        case 'T': return 4;
        case 'U': return 4;
        default: mexErrMsgIdAndTxt("SeqTool_nt2int:input", "Unknown character detected. Only NXACGTU are allowed. No lower case.");
    }
}

mxChar int2nt(int nt) {
    switch (nt) {
        case 0: return 'N';
        case 1: return 'A';
        case 2: return 'C';
        case 3: return 'G';
        case 4: return 'T';
        default: mexErrMsgIdAndTxt("SeqTool_int2nt:input", "Unknown int detected. Only integers 0-4 are allowed. 0-N, 1-A, 2-C, 3-G, 4-T.");
    }
}

void fixSeqDNA(char *pSeq) {
    for (size_t j = 0; j < std::char_traits<char>::length(pSeq); j++) {
        switch (pSeq[j]) {
            case 'A': break;
            case 'C': break;
            case 'G': break;
            case 'T': break;
            case 'N': break;
            case 'a': pSeq[j] = 'A'; break;
            case 'c': pSeq[j] = 'C'; break;
            case 'g': pSeq[j] = 'G'; break;
            case 't': pSeq[j] = 'T'; break;
            case 'u': pSeq[j] = 'T'; break;
            case 'U': pSeq[j] = 'T'; break;
            default: pSeq[j] = 'N'; break;
        }
    }
}

void fixSeqRNA(char *pSeq) {
    for (size_t j = 0; j < std::char_traits<char>::length(pSeq); j++) {
        switch (pSeq[j]) {
            case 'A': break;
            case 'C': break;
            case 'G': break;
            case 'U': break;
            case 'N': break;
            case 'a': pSeq[j] = 'A'; break;
            case 'c': pSeq[j] = 'C'; break;
            case 'g': pSeq[j] = 'G'; break;
            case 'u': pSeq[j] = 'U'; break;
            case 't': pSeq[j] = 'U'; break;
            case 'T': pSeq[j] = 'U'; break;
            default: pSeq[j] = 'N'; break;
        }
    }
}