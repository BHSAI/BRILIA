/*  
calcAlignScoreMEX takes a logical MxN matrix returned from a sequence
comparison resuls (SeqA == SeqB) and calculates the alignment score per
each row using the following equation:
SUM((consecutive matches)^2) - SUM((consecutive mismatches)^2) 

  Score = calcAlignScoreMEX(MatchResults)

  Score = calcAlignScoreMEX(MatchResults, AllowedMiss)

  Score = calcAlignScoreMEX(MatchResults, AllowedMiss, PenaltySide)
 
  INPUT
    MatchResults: 1xN logical array
    AllowedMiss: number of single-0 gap (eg, 1 0 1) that can be ignored 
      from left to right when calculating the consecutive match segment.
      EX: 1 0 1, AllowedMiss = 1 will return Score = 2^2 = 4.
      EX: 1 0 1, AllowedMiss = 0 will return Score = 1 - 1 + 1 = 1;
    PenaltySide: 'n', 'r', 'l', 'b' subtract penalty^2 for mismatched edges

  OUTPUT
    Score: the alignment score

  EXAMPLE
    MatchResults = [1 0 1 0 1 0 1 1 1 1]>0;
    AllowedMiss = 3;
    Score = calcAlignScoreMEX(MatchResults, AllowedMiss);
          =  49;
    Score = calcAlignScoreMEX(MatchResults, 0);
          =  16;
 */

#include "AlignTool.hpp"

void mexFunction(int nlhs,        mxArray *plhs[],
                 int nrhs, const  mxArray *prhs[]) {

    if (nrhs < 1 || nrhs > 3) {
        mexErrMsgIdAndTxt("calcAlignScoreMEX:nrhs", "Incorrect number of inputs. Min is 1. Max is 3.");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("calcAlignScoreMEX:nlhs", "Too many outputs. Max is 1.");
    }
    if (!mxIsLogical(prhs[0]) || mxGetM(prhs[0]) > 1) {
        mexErrMsgIdAndTxt("calcAlignScoreMEX:prhs", "Input1: MatchResults must be a 1xN logical array.");
    }
    if (nrhs >= 2 && !mxIsDouble(prhs[1])) {
        mexErrMsgIdAndTxt("calcAlignScoreMEX:prhs", "Input2: AllowedMiss must be an integer >= 0.");
    }
    if (nrhs >= 3 && !mxIsChar(prhs[2])) {
        mexErrMsgIdAndTxt("calcAlignScoreMEX:prhs", "Input3: PenaltySide must be 'n' (none), 'l' (left), 'r' (right), 'b' (both).");
    }
    
    bool *pMatch = mxGetLogicals(prhs[0]);
    mwSize Len = mxGetN(prhs[0]);   
    double AllowedMiss = nrhs >= 2 ? mxGetScalar(prhs[1]) : 0;
    mxChar PenaltySide = nrhs >= 3 ? *mxGetChars(prhs[2]) : 'n';
    
    plhs[0] = mxCreateDoubleScalar(calcAlignScore(pMatch, Len, AllowedMiss, PenaltySide));
}