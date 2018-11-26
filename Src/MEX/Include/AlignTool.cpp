/* Variable Descriptions
 *   ----------------------------------------------------------------------------------------------------------
 *   MissRate: Fraction of 1-gap-miss allowed for score calculation
 *   Alphabet: 'n', 'a', 'r' for nucleotide, amino acid, or random char matching
 *   ExactMatch: 'y', 'n' for exact match (same as cmprSeq) or alignment match
 *   TrimSide: 'l', 'r' for triming matches until 3/4 are found
 *   PenaltySide: 'l' (left), 'r' (right), or 'b' (both) side "leading" consec. misses reduce score by M^2.
 *   PreferSide: 'l', r', 'm' for tie-breaking alignment score to get left-match, right-match, or middle-match
 *   AllowedMiss: # of 1-gap-miss (101 but not 1001) to NOT be counted as a miss. 
 *  
 * NOTE: speed is faster with lower case options, as these uses if statements to check for parameter.
 */

#include "AlignTool.hpp"
#include <ctype.h>
#include <math.h>
#include <limits>

// Find the first occurence of a 1 (0-based index)
int findFirstMatch(bool *pMatch, mwSize Len) {
    for (int s = 0; s < Len; s++) {
        if (pMatch[s]) { return s; }
    }
    return -1; // no match
}

// Find the last occurence of a 1 (0-based index)
int findLastMatch(bool *pMatch, mwSize Len) {
    for (int e = Len-1; e >= 0; e--) {  
        if (pMatch[e]) { return e; }
    }
    return -1; // no match
}

// Compute the alignment score from a bool[] alignment result. 
// Score = sum(ConsecHits^2)-sum(ConsecMiss)^2-LeftMiss^2-RightMiss^2
double calcAlignScore(bool *pMatch, mwSize Len, double AllowedMiss, mxChar PenaltySide) {
    double Score = 0, Hits = 0, Miss = 0, PenaltyLHS = 0, PenaltyRHS = 0; // MUST initialize 0, or bugs occur (compiler bug)   
    int s = findFirstMatch(pMatch, Len);
    if (s < 0) { return - (double) (Len*Len); } //no need to check e<0 if s<0.
    int e = findLastMatch(pMatch, Len);
    
    for (int i = s; i <= e; i++) {
        if (pMatch[i]) {
            Hits++;
            if (Miss > 0) {
                Score -= Miss*Miss;
                Miss = 0;
            }
        } else {
            if (AllowedMiss > 0 && Miss == 0 && i+1 < Len && pMatch[i+1]) {
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
    
    switch (tolower(PenaltySide)) {
        case 'l':
            PenaltyLHS = s*s;
            break;
        case 'r':
            PenaltyRHS = (Len-e)*(Len-e);
            break;
        case 'b':
            PenaltyLHS = s*s;
            PenaltyRHS = (Len-e)*(Len-e);            
            break;
    }
    return (Score + Hits*Hits - Miss*Miss - PenaltyLHS - PenaltyRHS);
}

// Align SeqA and SeqB, returning an alignment information structure.
void alignSeq(mxChar *pSeqA, mxChar *pSeqB, mwSize LenA, mwSize LenB, double MissRate, mxChar Alphabet, mxChar ExactMatch, mxChar TrimSide, mxChar PenaltySide, mxChar PreferSide, align_info &AI) {
    bool pMatch[LenB > LenA ? LenB : LenA]; // only has to be greater of 2 sequence length. recycle array for match tracking.
    alignSeq(pSeqA, pSeqB, LenA, LenB, MissRate, Alphabet, ExactMatch, TrimSide, PenaltySide, PreferSide, AI, pMatch);
}
void alignSeq(mxChar *pSeqA, mxChar *pSeqB, mwSize LenA, mwSize LenB, double MissRate, mxChar Alphabet, mxChar ExactMatch, mxChar TrimSide, mxChar PenaltySide, mxChar PreferSide, align_info &AI, bool *pMatch) {
    if (LenA < 1 || LenB < 1) { return; } //Nothing to align 
    double AllowedMiss = 0;
    if (ExactMatch == 'y' || ExactMatch == 'Y') {
        mwSize Len = LenB > LenA ? LenA : LenB;
        cmprSeq(pSeqA, pSeqB, Len, 'n', pMatch);
        trimMatchResults(pMatch, Len, TrimSide);
        AllowedMiss = round(MissRate * Len);

        AI.Score  = calcAlignScore(pMatch, Len, AllowedMiss, PenaltySide);
        AI.BShift = 0;  //Reset in case some other function is recycling align_info
        AI.Match  = 0;  
        AI.MatchS = findFirstMatch(pMatch, Len);
        if (AI.MatchS > -1) { 
            AI.MatchE = findLastMatch(pMatch, Len);
            for (int z = AI.MatchS; z <= AI.MatchE; z++) { 
                if ( pMatch[z] ) { AI.Match++; };
            }
        } else {
            AI.MatchE = -1; //No match default
            AI.Score = - (double) Len*Len; //Maximum penalty allowed
        }

    } else { //ExactMatch == 'n';
        mxChar *pSeqL, *pSeqS; //L for longer seq, S for shorter seq
        mwSize LenL = 0, LenS = 0;
        if (LenB > LenA) {
            pSeqL = pSeqB;
            pSeqS = pSeqA;
            LenL = LenB;
            LenS = LenA;
        } else {
            pSeqL = pSeqA;
            pSeqS = pSeqB;
            LenL = LenA;
            LenS = LenB;
        }
        
        double AllowedMiss = 0;  // for determing allowed nt misses in alignment
        mwSize Len = 0;          // overlapping length of SeqL & SeqS
        int P = LenL + LenS - 2; // pos for pScore array
        double pScore[LenL + LenS - 1];
        
        //Reset align_info in case another program is recylcing AI 
        AI.Score = -1; //start at -1 because the first, worst you can do is Len = 1, all miss, so Score = - Len^2 = -1. 
        AI.Match =  0;
        AI.BShift = 0;
        AI.MatchS = 0;
        AI.MatchE = 0;
        
        //SeqS 1st nt is over the right edge of SeqL
        int L = LenL-1; //pos for SeqL
        int S = 0;      //pos for SeqS
        for (Len = 1; Len < LenS; Len++, L--, P--) {
            AllowedMiss = round(MissRate * Len);
            cmprSeq(&pSeqL[L], &pSeqS[S], Len, 'n', pMatch);
            trimMatchResults(pMatch, Len, TrimSide);
            pScore[P] = calcAlignScore(pMatch, Len, AllowedMiss, PenaltySide);
            if (pScore[P] > AI.Score) { AI.Score = pScore[P]; } //only AI.Score needs to be update to get highest score
        }

        //SeqS is over SeqL. Len should be same as LenS.
        for (L; L >= 0; L--, P--) {
            AllowedMiss = round(MissRate * Len);
            cmprSeq(&pSeqL[L], &pSeqS[S], Len, 'n', pMatch);
            trimMatchResults(pMatch, Len, TrimSide);
            pScore[P] = calcAlignScore(pMatch, Len, AllowedMiss, PenaltySide);
            if (pScore[P] > AI.Score) { AI.Score = pScore[P]; }
        }
        
        //SeqS lst nt is over the left edge of SeqL
        L = 0;
        S = 1;
        for (Len = LenS-1; Len >= 1; Len--, S++, P--) {
            AllowedMiss = round(MissRate * Len);
            cmprSeq(&pSeqL[L], &pSeqS[S], Len, 'n', pMatch);
            trimMatchResults(pMatch, Len, TrimSide);
            pScore[P] = calcAlignScore(pMatch, Len, AllowedMiss, PenaltySide);
            if (pScore[P] > AI.Score) { AI.Score = pScore[P]; }
        }

        //Find max score location
        if (PreferSide == 'l' || PreferSide == 'L') { //Find left highest
            for (P = 0; P <= LenS+LenL-2; P++) {
                if (fabs(pScore[P] - AI.Score) < std::numeric_limits<double>::epsilon()) { break; }
            }
        } else if (PreferSide == 'r' || PreferSide == 'R') { //Find right highest
            for (P = LenS+LenL-2; P >= 0; P--) {
                if (fabs(pScore[P] - AI.Score) < std::numeric_limits<double>::epsilon()) { break; }
            }
        } else { //Find middle highest
            int Loc[LenS+LenL-1];
            int q = 0;
            for (P = 0; P <= LenS+LenL-2; P++) {
                if (fabs(pScore[P] - AI.Score) < std::numeric_limits<double>::epsilon()) { 
                    Loc[q++] = P; 
                }
            }
            P = Loc[(int) round(q/2)];
        }
   
        //Determine the AI parameters
        int BShift = (P - (int) LenS + 1);
        AI.BShift = LenB > LenA ? -BShift : BShift;

        if (AI.BShift >= 0) {
            Len = LenA - AI.BShift < LenB ? LenA - AI.BShift : LenB;
            cmprSeq(&pSeqA[AI.BShift],  &pSeqB[0], Len, Alphabet, pMatch);
        } else { 
            Len = LenA - AI.BShift < LenB ? LenA : LenB + AI.BShift;
            cmprSeq(&pSeqA[0], &pSeqB[-AI.BShift], Len, Alphabet, pMatch);
        }
        trimMatchResults(pMatch, Len, TrimSide);

        AI.Match  = 0;
        AI.MatchS = findFirstMatch(pMatch, Len);
        if (AI.MatchS > -1) { 
            AI.MatchE = findLastMatch(pMatch, Len);
            for (int z = AI.MatchS; z <= AI.MatchE; z++) { 
                if ( pMatch[z] ) { AI.Match++; };
            }
            AI.MatchS += fabs(AI.BShift);
            AI.MatchE += fabs(AI.BShift);
        } else {
            AI.Score = - (double) Len*Len; //Maximum penalty allowed
            AI.MatchS += fabs(AI.BShift);
            AI.MatchE = -1; //No match default
        }
    }
} 

// Compares SeqA and SeqB and updates a boolean vector of match/miss.
void cmprSeq(mxChar *pSeqA, mxChar *pSeqB, mwSize Len, mxChar Alphabet, bool *pMatch) { // overloaded: compare SeqA and SeqB WITHOUT creating a bool *palignment result
    switch (tolower(Alphabet)) {
        case 'n':
            for (mwSize i = 0; i < Len; i++) {
                pMatch[i] = (pSeqA[i] == pSeqB[i] || pSeqA[i] == 'N' || pSeqB[i] == 'N');
            }
            break;
        case 'a':
            for (mwSize i = 0; i < Len; i++) {
                pMatch[i] = (pSeqA[i] == pSeqB[i] || pSeqA[i] == 'X' || pSeqB[i] == 'X');
            }
            break;
        default:
            for (mwSize i = 0; i < Len; i++) {
                pMatch[i] = (pSeqA[i] == pSeqB[i]);
            }        
            break;
    }
}

// Trims (set to false) bool[] from the left and/or right side until 3/4 matches are found.
void trimMatchResults(bool *pMatch, mwSize Len, mxChar TrimSide) {
    TrimSide = tolower(TrimSide);
    if (Len <= 4 || TrimSide == 'n') { return; }
    
    if (TrimSide == 'l' || TrimSide == 'b') {
        int Sum = 0, i = 0;
        for (i; i < 4; i++) { //get initial 1st 4 sum
            if (pMatch[i]) { Sum++; } 
        }
        for (i; i < Len; i++) { //i is 4 and going up
            if (Sum >= 3) { break; }
            if (pMatch[i]) { Sum++; }
            if (pMatch[i-4]) { 
                Sum--;
                pMatch[i-4] = false;
            }
        }
        if (i == Len && Sum < 3) { //set all to 0, nothing is good
            for (i = i-4; i < Len; i++) {
                pMatch[i] = false;
            }
        }
    }

    if (TrimSide == 'r' || TrimSide == 'b') {
        int Sum = 0, i = Len;
        for (i; i > Len-4; i--) {//Initial 1st 4 sum 
            if (pMatch[i-1]) { Sum++; }
        }
        for (i; i >= 1; i--) {
            if (Sum >= 3) { break; }
            if (pMatch[i-1]) { Sum++; } 
            if (pMatch[i+3]) { 
                Sum--; 
                pMatch[i+3] = false; 
            }
        }
        if (i == 0 && Sum < 3) { //set all to 0, nothing is good
            for (i = 3; i >= 0; i--) {
                pMatch[i] = false;
            }
        }
    }
}

mxArray *buildAlignment(mxChar *pSeqA, mxChar *pSeqB, mwSize LenA, mwSize LenB, align_info AI) {
    mwSize Dims[2] = {3, 0};
    int a0 = 0, b0 = 0;
    if (AI.BShift >= 0) {
        Dims[1] = LenA >= LenB + AI.BShift ? LenA : LenB + AI.BShift;
        b0 =  AI.BShift;
    } else {
        Dims[1] = LenA - AI.BShift >= LenB ? LenA - AI.BShift : LenB;
        a0 = -AI.BShift;
    }

    mxArray *pAlign = mxCreateCharArray(2, Dims);   
    mxChar *pAlignStr = mxGetChars(pAlign);
    mwSize Len = AI.MatchS < 0 ? 0 : AI.MatchE - AI.MatchS + 1;

    for (mwSize a = 0; a < LenA; a++) {
        pAlignStr[0 + 3*(a+a0)] = pSeqA[a];
    }
    for (mwSize b = 0; b < LenB; b++) {
        pAlignStr[2 + 3*(b+b0)] = pSeqB[b];
    }
    
    if (AI.Match > 0) {
        bool pMatch[Len];
        cmprSeq(&pSeqA[AI.MatchS - a0], &pSeqB[AI.MatchS - b0], Len, 'n', pMatch);
        int t = findFirstMatch(pMatch, Len);
        for (int m = AI.MatchS; m <= AI.MatchE; m++, t++) {
            pAlignStr[1 + 3*m] = pMatch[t] ? '|' : ' ';
        }
    }
    return pAlign;
}

int countMatch(bool *pMatch, mwSize Len) {
    int Ct;
    for (int k = 0; k < Len; k++) {
        if (pMatch[k]) { Ct++; }
    }
    return Ct;
}