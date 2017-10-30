#include "mex.h"
#include <math.h>
#include <string>
/*
SeqA = 'ACGTACGT';
SeqB = 'ACGT';       
MissRate = 0;
Alphabet = 'n';
alignSeqMEX(SeqA, SeqB, MissRate, Alphabet);
*/

/* cmprSeq will compare SeqA to SeqB, and adjust a bool[0:Len-1] to 1 for 
 * match, 0 otherwise.
 */
void cmprSeq(mxChar *pSeqA, mxChar *pSeqB, mwSize Len, mxChar Alphabet, bool *pMatch) {
    switch (tolower(Alphabet)) {
        case 'n':
            for (mwSize i = 0; i < Len; i++) {
                if (pSeqA[i] == pSeqB[i] ||
                    pSeqA[i] == 'N' || pSeqA[i] == 'X' || 
                    pSeqB[i] == 'N' || pSeqB[i] == 'X') {
                    pMatch[i] = true;
                } else {
                    pMatch[i] = false;
                }
            }
            break;
        case 'a':
            for (mwSize i = 0; i < Len; i++) {
                if (pSeqA[i] == pSeqB[i] || pSeqA[i] == 'X' || pSeqB[i] == 'X') {
                    pMatch[i] = true;
                } else {
                    pMatch[i] = false;
                }
            }
            break;
        default:
            for (mwSize i = 0; i < Len; i++) {
                if (pSeqA[i] == pSeqB[i]) {
                    pMatch[i] = true;
                } else {
                    pMatch[i] = false;
                }
            }        
            break;
    }
} //cmprSeq

/* trimMatchResults will trim bool[0:Len-1] on the Left and Right side  
 * until a 3 out of 4 match zone is found.
 */
void trimMatchResults(bool *pMatch, mwSize Len, mxChar TrimSide) {
    TrimSide = tolower(TrimSide);
    if (Len <= 4 || TrimSide == 'n') {
        return;
    }
    
    int Sum;
    mwSize i;

    if (TrimSide == 'l' || TrimSide == 'b') {
        Sum = 0;
        for (i = 0; i < 4; i++) {//Initial 1st 4 sum
            if (pMatch[i]) { 
                Sum++; 
            }
        }
        for (i; i < Len; i++) {
            if (Sum >= 3) { 
                break; 
            }
            if (pMatch[i-4]) { 
                Sum--; 
                pMatch[i-4] = 0; 
            } 
            if (pMatch[i]) {
                Sum++;
            }
        }
    }

    if (TrimSide == 'r' || TrimSide == 'b') {
        Sum = 0;
        for (i = Len; i > Len-4; i--) {//Initial 1st 4 sum (NOTE: mwSize is unsigned, so i will NEVER be negative.)
            if (pMatch[i-1]) {
                Sum++;
            }
        }
        for (i; i >= 1; i--) {
            if (Sum >= 3) {
                break;
            }
            if (pMatch[i+3]) {
                Sum--;
                pMatch[i+3] = 0;
            }
            if (pMatch[i-1]) {
                Sum++;
            }
        }
    }
} //trimMatchResults


/* calcAlignScore will computer the alignment score of a bool[0:Len-1] 
 * using the sum(ConsecHits^2) - sum(ConsecMiss)^2. AllowedMiss allows for 
 * a one-gap-miss (101 but not 1001) to be treated as a 0 Hit. PenaltySide 
 * will subtract Miss^2 of leading ('l') or trailing ('r') 0's.
 */
void calcAlignScore(bool *pMatch, mwSize Len, int AllowedMiss, mxChar TrimSide, mxChar PenaltySide, int *pScore) {
    int Score = 0;
    int Hits  = 0;
    int Miss  = 0;
    int LeftPenalty  = 0;
    int RightPenalty = 0;
    mwSize s, e, i;
    
    trimMatchResults(pMatch, Len, TrimSide);
    
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
            if (i+1 < Len && pMatch[i+1] && AllowedMiss > 0) {
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
    *pScore = Score + Hits*Hits - Miss*Miss - LeftPenalty - RightPenalty;
} //calcAlignScore  

/* slignSeqMEX is the main function for doing alignments between SeqA and 
 * SeqB based on the other parameters given.
 */
void alignSeqMEX(mxChar *pSeqA, mxChar *pSeqB, mwSize LenA, mwSize LenB, double MissRate, mxChar Alphabet, mxChar ExactMatch, mxChar TrimSide, mxChar PreferSide, mxChar PenaltySide, int *pMaxScore, int *pBshift) {
    if (LenA == 0 || LenB == 0) {
        return;
    }
            
    if (tolower(ExactMatch) == 'y') {
        mwSize MinLen = LenA < LenB ? LenA : LenB;
        bool pMatch[MinLen];
        cmprSeq(pSeqA, pSeqB, MinLen, 'n', pMatch);
        int AllowedMiss = (int) round(MissRate * MinLen);
        calcAlignScore(pMatch, MinLen, AllowedMiss, TrimSide, PenaltySide, pMaxScore);
        *pBshift = 0;
        
    } else { //ExactMatch == 'n';
        //Want SeqA to be the longer one, SeqB is the shorter one
        bool IsAgtB = LenA > LenB ? true : false;
        if (!IsAgtB) {
            mxChar *pTemp = pSeqA;
            pSeqA = pSeqB;
            pSeqB = pTemp;
            mwSize TempLen = LenA;
            LenA = LenB;
            LenB = TempLen;
        }

        mwSize Len, a, b, q = LenA+LenB-2; // pos counter for SeqA, SeqB, and pScore
        int  AllowedMiss = 0;
        bool pMatch[LenB];
        int  pScore[LenA+LenB-1];
        
        //SeqB 1st nt is over the right edge of SeqA
        a = LenA-1;
        b = 0;
        for (Len = 1; Len < LenB; Len++, a--, q--) {
            AllowedMiss = (int) round(MissRate * Len);
            cmprSeq(&pSeqA[a], &pSeqB[b], Len, 'n', pMatch);
            calcAlignScore(pMatch, Len, AllowedMiss, TrimSide, PenaltySide, &pScore[q]);
            if (pScore[q] > *pMaxScore) {
                *pMaxScore = pScore[q];
            }
        }
        
        //SeqB is over SeqA
        for (a; a+1 >= 1; a--, q--) {
            cmprSeq(&pSeqA[a], &pSeqB[b], LenB, 'n', pMatch);
            calcAlignScore(pMatch, LenB, AllowedMiss, TrimSide, PenaltySide, &pScore[q]);
            if (pScore[q] > *pMaxScore) {
                *pMaxScore = pScore[q];
            }
        }

        //SeqB lst nt is over the left edge of SeqA
        a = 0;
        b = 1;
        for (Len = LenB-1; Len >= 1; Len--, b++, q--) {
            AllowedMiss = (int) round(MissRate * Len);
            cmprSeq(&pSeqA[a], &pSeqB[b], Len, 'n', pMatch);
            calcAlignScore(pMatch, Len, AllowedMiss, TrimSide, PenaltySide, &pScore[q]);
            if (pScore[q] > *pMaxScore) {
                *pMaxScore = pScore[q];
            }
        }

        //Find max score location
        switch (tolower(PreferSide)) {
            case 'l': //Find left highest
                for (q = 0; q < LenB+LenA-1; q++) {
                    if (pScore[q] == *pMaxScore) {
                        break;
                    }
                }
                break;
            case 'r': //Find right highest
                for (q = LenB+LenA-1; q-- > 0;) {
                    if (pScore[q] == *pMaxScore) {
                        break;
                    }
                }
                break;
            default:  //Find middle highest
                int Loc[LenB+LenA-1];
                int p = 0;
                for (q = 0; q < LenB+LenA-1; q++) {
                    if (pScore[q] == *pMaxScore) 
                        Loc[p++] = q;
                }
                q = Loc[(int) floor(p/2)];
                break;
        }

        //Outputs flipped depending on if LenA > LenB
        if (IsAgtB) {
            *pBshift =   q - (int) LenB + 1; 
        } else {
            *pBshift = -(q - (int) LenB + 1);  
        }
    }
} 

//Following functions are for assembling the outputs
void getAlignmentInfo(bool *pMatch, mxChar *pSeqA, mxChar *pSeqB, mwSize LenA, mwSize LenB, mxChar Alphabet, mxChar TrimSide, int Bshift, int *pScore, int *pStartAt, int *pMatchAt) {
    mwSize Len;
    if (Bshift >= 0) {
        Len = LenA - Bshift < LenB ? LenA - Bshift : LenB;
        cmprSeq(&pSeqA[Bshift],  &pSeqB[0], Len, Alphabet, pMatch);
    } else { 
        Len = LenA - Bshift < LenB ? LenA : LenB + Bshift;
        cmprSeq(&pSeqA[0], &pSeqB[-Bshift], Len, Alphabet, pMatch);
    }
    
    trimMatchResults(pMatch, Len, TrimSide);

    int s, e;
    bool HasMatch = false;
    for (s = 0; s < Len; s++) {
        if (pMatch[s]) {
            HasMatch = true;
            break;
        }
    }
    if (!HasMatch) {
        return;
    }
    for (e = Len; e >= 1; e--) {
        if (pMatch[e-1]) {
            break;
        }
    }
    
    for (int i = s; i < e; i++) {
        if (pMatch[i]) {
            pScore[0]++;
        }
    }

    pStartAt[0] = 1; //Matlab doesn't use the 0-index, so 1.
    pStartAt[1] = Bshift >= 0 ? Bshift + 1 : Bshift; //Note: StartAt = 1 is s = 0, but StartAt = -1 is s = -1

    pMatchAt[0] = s + abs(Bshift) + 1; //Match of entire alignment text, add 1 for matlab index convention.
    pMatchAt[1] = e + abs(Bshift); //Note: same as "1 + e + abs(Bshift) - 1", which adjusts for matlab conv. and inclusive end.
}
        
void mexFunction(int nlhs,        mxArray *plhs[],
                 int nrhs, const  mxArray *prhs[]) {

    if (nrhs < 2)
        mexErrMsgIdAndTxt("alignSeqMEX:input", "Not enough inputs.");

    if (nlhs > 4)
        mexErrMsgIdAndTxt("alignSeqMEX:output", "Too many outputs.");
    
    if (!mxIsChar(prhs[0]) || mxGetM(prhs[0]) > 1)
        mexErrMsgIdAndTxt("alignSeqMEX:input", "Input1: SeqA, must be a 1xN char array.");
    
    if (nrhs >= 2 && (!mxIsChar(prhs[1]) || mxGetM(prhs[1]) > 1))
        mexErrMsgIdAndTxt("alignSeqMEX:input", "Input2: SeqB, must be a 1xN char array.");
    
    if (nrhs >= 3 && !mxIsDouble(prhs[2]))
        mexErrMsgIdAndTxt("alignSeqMEX:input", "Input3: MissRate, must be a scalar between 0.0 to 1.0.");
    
    if (nrhs >= 4 && !mxIsChar(prhs[3]))
        mexErrMsgIdAndTxt("alignSeqMEX:input", "Input4: Alphabet, must be a char 'n' (nucleotide), 'a' (amino acid), or other (basic char).");
    
    if (nrhs >= 5 && !mxIsChar(prhs[4]))
        mexErrMsgIdAndTxt("alignSeqMEX:input", "Input5: ExactMatch, must be a char 'n' or 'y'.");

    if (nrhs >= 6 && !mxIsChar(prhs[5]))
        mexErrMsgIdAndTxt("alignSeqMEX:input", "Input6: TrimSide, must be a char 'n' (none), 'l' (left), 'r' (right), 'b' (both).");
    
    if (nrhs >= 7 && !mxIsChar(prhs[6]))
        mexErrMsgIdAndTxt("alignSeqMEX:input", "Input7: PenaltySide, must be a char 'n' (none), 'l' (left), 'r' (right), 'b' (both).");

    if (nrhs >= 8 && !mxIsChar(prhs[7]))
        mexErrMsgIdAndTxt("alignSeqMEX:input", "Input8: PreferSide, must be a char 'n' (none), 'l' (left), 'r' (right).");

    mxChar *pSeqA, *pSeqB;
    mwSize LenA, LenB; 
    mxChar Alphabet, ExactMatch, TrimSide, PenaltySide, PreferSide;
    double MissRate;

    pSeqA = mxGetChars(prhs[0]);
    pSeqB = mxGetChars(prhs[1]);
    LenA  = mxGetN(prhs[0]);
    LenB  = mxGetN(prhs[1]);

    if (nrhs >= 3) {
        MissRate = mxGetScalar(prhs[2]);
    } else {
        MissRate = 0;
    }

    if (nrhs >= 4) {
        Alphabet = *mxGetChars(prhs[3]);
    } else {
        Alphabet = 'n';
    }

    if (nrhs >= 5) {
        ExactMatch = *mxGetChars(prhs[4]);
    } else {
        ExactMatch = 'n';
    }

    if (nrhs >= 6) {
        TrimSide = *mxGetChars(prhs[5]);
    } else {
        TrimSide = 'n';
    }

    if (nrhs >= 7) {
        PenaltySide = *mxGetChars(prhs[6]);
    } else {
        PenaltySide = 'n';
    }

    if (nrhs >= 8) {
        PreferSide = *mxGetChars(prhs[7]);
    } else {
        PreferSide = 'n';
    }
    
    int Bshift = 0;
    int pScore[2] = {0, INT_MIN}; //[MatchCt; MaxScore]
    int pStartAt[2] = {0};
    int pMatchAt[2] = {0}; 
    bool pMatch[LenA+LenB];
    
    alignSeqMEX(pSeqA, pSeqB, LenA, LenB, MissRate, Alphabet, ExactMatch, TrimSide, PreferSide, PenaltySide, &pScore[1], &Bshift); 
    getAlignmentInfo(pMatch, pSeqA, pSeqB, LenA, LenB, Alphabet, TrimSide, Bshift, pScore, pStartAt, pMatchAt);

    if (nlhs >= 1) {
        plhs[0] = mxCreateDoubleMatrix(2, 1, mxREAL);
        double *pScoreOut = mxGetPr(plhs[0]);
        for (int k = 0; k < 2; k++) {
            pScoreOut[k] = pScore[k];
        }
    }

    if (nlhs >= 2) {
        plhs[1] = mxCreateDoubleMatrix(2, 1, mxREAL);
        double *pStartAtOut = mxGetPr(plhs[1]);
        for (int k = 0; k < 2; k++) {
            pStartAtOut[k] = pStartAt[k];
        }
    }

    if (nlhs >= 3) {
        plhs[2] = mxCreateDoubleMatrix(2, 1, mxREAL);
        double *pMatchAtOut = mxGetPr(plhs[2]);
        for (int k = 0; k < 2; k++) {
            pMatchAtOut[k] = pMatchAt[k];
        }
    }

    if (nlhs >= 4) { //Make the alignment
        mwSize Dims[2] = {3, 1};
        if (Bshift >= 0) {
            Dims[1] = LenA >= LenB + Bshift ? LenA : LenB + Bshift;
            plhs[3] = mxCreateCharArray(2, Dims);
            mxChar *pAlignStr = mxGetChars(plhs[3]);
            for (mwSize a = 0; a < LenA; a++) {
                pAlignStr[0+3*a] = pSeqA[a];
            }
            for (mwSize b = 0; b < LenB; b++) {
                pAlignStr[2+3*(b+Bshift)] = pSeqB[b];
            }
            mwSize t = 0;
            mwSize Len = LenA - Bshift < LenB ? LenA - Bshift : LenB;
            for (mwSize m = Bshift; m < Bshift + Len; m++, t++) {
                pAlignStr[1+3*m] = pMatch[t] ? '|' : ' ';
            }
        } else {
            Dims[1] = LenA - Bshift >= LenB ? LenA - Bshift : LenB;
            plhs[3] = mxCreateCharArray(2, Dims);
            mxChar *pAlignStr = mxGetChars(plhs[3]);
            for (mwSize a = 0; a < LenA; a++) {
                pAlignStr[0+3*(a-Bshift)] = pSeqA[a];
            }
            for (mwSize b = 0; b < LenB; b++) {
                pAlignStr[2+3*b] = pSeqB[b];
            }
            mwSize t = 0;
            mwSize Len = LenA - Bshift < LenB ? LenA : LenB + Bshift;
            for (mwSize m = -Bshift; m < -Bshift + Len; m++, t++) {
                pAlignStr[1+3*m] = pMatch[t] ? '|' : ' ';
            }
        }
    }
}