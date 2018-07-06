/*  
cmprSeqMEX will compare two sequences and return a logical array where
1 = match and 0 = mismatch. 
  
  Match = cmprSeqMEX(SeqA, SeqB, Alphabet)
   
  INPUT
    SeqA: 1st sequence (CASE SENSITIVE!)
    SeqB: 2nd sequence (CASE SENSITIVE!)
    Alphabet ['n', 'a', 'r']: nucleotide, amino acid, or random sequenc
      'n': X and N are wildcard for DNA/RNA (default)
      'a': X is wildcard for AA
      'r': no wildcards for character
    
  OUTPUT
    Match: 1-row logical array of matches (1) and mismatch(0)

  EXAMPLE
    SeqA = 'ACGTXXNNACGT';
    SeqB = 'ACGTACGTACGT';

    Match = cmprSeqMEX(SeqA, SeqB, 'n')
        =  1  1  1  1  1  1  1  1  1  1  1  1
    Match = cmprSeqMEX(SeqA, SeqB, 'a')
        =  1  1  1  1  1  1  1  1  1  1  1  1
    Match = cmprSeqMEX(SeqA, SeqB, 'r')
        =  1  1  1  1  0  0  0  0  1  1  1  1
*/

#include "mex.h"

void cmprSeqNT(const mxChar *CharA, const mxChar *CharB, mwSize M, bool *Match) {
    mwSize i;
    for (i = 0; i < M; i++) {
        if (CharA[i] == CharB[i] ||
            CharA[i] == 'N' || CharA[i] == 'X' || 
            CharB[i] == 'N' || CharB[i] == 'X') {
 
            Match[i] = true;
        }
    }
}

void cmprSeqAA(const mxChar *CharA, const mxChar *CharB, mwSize M, bool *Match) {
    mwSize i;
    for (i = 0; i < M; i++) {
        if (CharA[i] == CharB[i] ||
            CharA[i] == 'X' || 
            CharB[i] == 'X') {
 
            Match[i] = true;
        }
    }
}

void cmprSeq(const mxChar *CharA, const mxChar *CharB, mwSize M, bool *Match) {
    mwSize i;
    for (i = 0; i < M; i++) {
        if (CharA[i] == CharB[i]) {
 
            Match[i] = true;
        }
    }
}

void mexFunction(int nlhs,        mxArray *plhs[],
                 int nrhs, const  mxArray *prhs[]) {

    if (nrhs < 2) {
        mexErrMsgIdAndTxt("cmprSeqMEX:nrhs", "Need >= 2 inputs.");
    }
    
    if (!mxIsChar(prhs[0])) { 
        mexErrMsgIdAndTxt("cmprSeqMEX:prhs", "Must be a char array.");
    }
    
    if (mxGetM(prhs[0]) > 1) {
        mexErrMsgIdAndTxt("cmprSeqMEX:prhs", "Must be a 1xN vector.");
    }
    
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("cmprSeqMEX:nlhs", "Need 1 output.");
    }

    mwSize NA = mxGetN(prhs[0]);
    mwSize NB = mxGetN(prhs[1]);
    if (NA != NB) {
        mexErrMsgIdAndTxt("cmprSeqMEX:prhs", "SeqA and SeqB must be same size.");
    }
    
    char Alphabet = 'n';
    if (nrhs >= 3) {
        Alphabet = *mxGetChars(prhs[2]);
    }
    
    plhs[0] = mxCreateLogicalMatrix(1, NA);
    bool *Match = mxGetLogicals(plhs[0]);

    mxChar *CharA = mxGetChars(prhs[0]);
    mxChar *CharB = mxGetChars(prhs[1]);
    switch (Alphabet) {
        case 'n': cmprSeqNT(CharA, CharB, NA, Match); break;
        case 'a': cmprSeqAA(CharA, CharB, NA, Match); break;
        case 'r': cmprSeq(CharA, CharB, NA, Match);   break;
        default : mexErrMsgIdAndTxt("cmprSeqMEX:prhs", "Unrecognized 2nd input. Must be 'n', 'a', or 'r'.");
    }                
}