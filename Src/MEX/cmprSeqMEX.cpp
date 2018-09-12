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
        =  1  1  1  1  1  1  0  0  1  1  1  1
    Match = cmprSeqMEX(SeqA, SeqB, 'r')
        =  1  1  1  1  0  0  0  0  1  1  1  1
*/

#include "AlignTool.hpp"

void mexFunction(int nlhs,        mxArray *plhs[],
                 int nrhs, const  mxArray *prhs[]) {

    if (nrhs < 2 || nrhs > 3) {
        mexErrMsgIdAndTxt("cmprSeqMEX:nrhs", "Incorrect number of inputs. Min is 2. Max is 3.");
    }
    
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("cmprSeqMEX:nlhs", "Too many outputs. Max is 1.");
    }

    if (!mxIsChar(prhs[0])  || !mxIsChar(prhs[1]) || 
        mxGetM(prhs[0]) > 1 || mxGetM(prhs[1]) > 1) { 
        mexErrMsgIdAndTxt("cmprSeqMEX:prhs", "First 2 inputs must be 1xN char arrays.");
    }
    
    mwSize LenA = mxGetN(prhs[0]);
    mwSize LenB = mxGetN(prhs[1]);
    mwSize Len = LenA >= LenB ? LenB : LenA;
    mxChar Alphabet = nrhs >= 3 ? *mxGetChars(prhs[2]) : 'n';
    
    plhs[0] = mxCreateLogicalMatrix(1, Len);
    bool *pMatch = mxGetLogicals(plhs[0]);
    mxChar *pSeqA = mxGetChars(prhs[0]);
    mxChar *pSeqB = mxGetChars(prhs[1]);
    
    cmprSeq(pSeqA, pSeqB, Len, Alphabet, pMatch);
}