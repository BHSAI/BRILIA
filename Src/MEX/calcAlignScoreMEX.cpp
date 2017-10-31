/*  
calcAlignScore takes a logical MxN matrix returned from a sequence
comparison resuls (SeqA == SeqB) and calculates the alignment score per
each row using the following equation:
SUM((consecutive matches)^2) - SUM((consecutive mismatches)^2) 

Score = calcAlignScore(MatchResults)

Score = calcAlignScore(MatchResults, AllowedMiss)

INPUT
  MatchResults: 1xN logical array
  AllowedMiss: number of single-0 gap (eg, 1 0 1) that can be ignored 
    from left to right when calculating the consecutive match segment.
    EX: 1 0 1, AllowedMiss = 1 will return Score = 2^2 = 4.
    EX: 1 0 1, AllowedMiss = 0 will return Score = 1 - 1 + 1 = 1;

OUTPUT
  Score: the alignment score

EXAMPLE
  MatchResults = [1 0 1 0 1 0 1 1 1 1];
  AllowedMiss = 3;
  Score = calcAlignScoreMEX(MatchResults, AllowedMiss);
    =  49;
  Score = calcAlignScoreMEX(MatchResults, 0);
    =  16;
 */

#include "mex.h"
#include <string>

/* calcAlignScore will computer the alignment score of a bool[0:Len-1] 
 * using the sum(ConsecHits^2) - sum(ConsecMiss)^2. AllowedMiss allows for 
 * a one-gap-miss (101 but not 1001) to be treated as a 0 Hit. PenaltySide 
 * will subtract Miss^2 of leading ('l') or trailing ('r') 0's.
 */
void calcAlignScore(bool *pMatch, mwSize Len, int AllowedMiss, mxChar PenaltySide, double *pScore) {
    int Score = 0;
    int Hits  = 0;
    int Miss  = 0;
    int LeftPenalty  = 0;
    int RightPenalty = 0;
    mwSize s, e, i;
    
    for (s = 0; s < Len; s++) {
        if (pMatch[s]) {
            break;
        }
    }
    for (e = Len; e >= 1; e--) {
        if (pMatch[e-1]) {
            break;
        }
    }
    for (i = s; i < e; i++) {
        if (pMatch[i]) {
            Hits++;
            if (Miss > 0) {
                Score -= Miss*Miss;
                Miss = 0;
            }
        } else {
            if (Miss == 0 && i+1 < Len && pMatch[i+1] && AllowedMiss > 0) {
                AllowedMiss--;
            } else {
                Miss++;
                if (Hits > 0) {
                    Score += Hits*Hits;
                    Hits = 0;
                }
            }
        }
    }
    
    PenaltySide = tolower(PenaltySide);
    if (PenaltySide == 'l' || PenaltySide == 'b') {
        LeftPenalty  = s*s;
    }
    if (PenaltySide == 'r' || PenaltySide == 'b') {
        RightPenalty = (Len - e)*(Len - e);
    }
    *pScore = (double) Score + Hits*Hits - Miss*Miss - LeftPenalty - RightPenalty;
} //calcAlignScore  

void mexFunction(int nlhs,        mxArray *plhs[],
                 int nrhs, const  mxArray *prhs[]) {

    if (nrhs < 1)
        mexErrMsgIdAndTxt("calcAlignScoreMEX:input", "Not enough inputs.");

    if (nlhs > 1)
        mexErrMsgIdAndTxt("calcAlignScoreMEX:output", "Too many outputs.");
    
    if (!mxIsLogical(prhs[0]) || mxGetM(prhs[0]) > 1)
        mexErrMsgIdAndTxt("calcAlignScoreMEX:input", "Input1: MatchResults must be a 1xN binary array.");
    
    if (nrhs >= 2 && !mxIsDouble(prhs[1]))
        mexErrMsgIdAndTxt("calcAlignScoreMEX:input", "Input2: AllowedMiss must be an integer >= 0.");

    if (nrhs >= 3 && !mxIsChar(prhs[2]))
        mexErrMsgIdAndTxt("calcAlignScoreMEX:input", "Input3: PenaltySide must be a char 'n' (none), 'l' (left), 'r' (right), 'b' (both).");

    bool *pMatch = mxGetLogicals(prhs[0]);
    mwSize N = mxGetN(prhs[0]);
    
    double AllowedMiss = 0;
    if (nrhs >= 2) {
        AllowedMiss = mxGetScalar(prhs[1]);
    }
    
    mxChar PenaltySide = 'n';
    if (nrhs >= 3) {
        PenaltySide = *mxGetChars(prhs[2]);
    }
    
    double Score = 0;
    calcAlignScore(pMatch, N, AllowedMiss, PenaltySide, &Score);
    plhs[0] = mxCreateDoubleScalar(Score);
}