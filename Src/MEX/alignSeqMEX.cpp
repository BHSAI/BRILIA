/*
alignSeqMEX will align sequence SeqA and SeqB using the BRILIA alignment
protocol. Consecutive matches add C^2 score and consecutive mismatches
subtract M^2 score, where C and M is the length of each contiguous
segment. If AllowedMiss is specified, then C segments can be elongated
only if there is a 1-letter miss, and if the number of AllowedMiss has not
been reached, going left to right. More information is provided in the
reference for BRILIA.

  Score = alignSeqMEX(SeqA, SeqB)

  [Score, StartAt] = alignSeqMEX(SeqA, SeqB)

  [Score, StartAt, MatchAt] = alignSeqMEX(SeqA, SeqB)

  [Score, StartAt, MatchAt, Alignment] = alignSeqMEX(SeqA, SeqB)

  ... = alignSeqMEX(SeqA, SeqB, MissRate, Alphabet, ExactMatch, TrimSide,
          PenaltySide, PreferSide)
  
  INPUTS
    SeqA: Character sequence. X = wildcard match. Z = do not match.
    SeqB: Character sequence. X = wildcard match. Z = do not match.
    MissRate [Integer >= 0]: How many missed nt/bp are allowed for
      elongating a contiguous C region for alignment score. Ex: MissRate
      0.02 = 2 of shorter seq length.
    Alphabet ['n' 'a' 'r']: Specify nt, aa, or random char. For nt, N and
      X char are wildcard. for aa, only X is wildcard. For random, no
      wildcards.
    ExactMatch ['n' 'y']: If SeqA and SeqB are same lengths and you
      want to do a simple, direct alignment, then use 'y'.
    TrimSide ['n' 'l','r' 'b']: Will perform a 3/4 match on the side:
      none, left, right, or both.
    PenaltySide ['n' 'l' 'r']: Will reduce alignment score if
      none, left, or right portion of SeqA overhangs SeqB and is unused.
      Penalty is (unused nt)^2.
    PreferSide ['n' 'l' 'r']: If a tie alignment is found, will favor the
      alignment closers to this side of SeqA.

  OUTPUTS
    Score: 2x1 matrix, Score(1) is the # of hits, and Score(2) is the
      aligment score of BRILIA.
    StartAt: 2x1 matrix, StartAt(1) is always 1 the SeqA start, and
      StartAt(2) is the SeqB start relative to the nth letter of SeqA. 
      EX: if StartAt(2) = 3, then SeqB aligns with A if its 1st letter
      start below the 3rd letter of SeqB. If StartAt(2) = -3, then the 1st
      letter of SeqB start left of SeqA's 1st letter by 3 letters (missing
      3 letters).
    MatchAt: 2x1 matrix showing the 1st and last location of the '|' in
      the Alignment char matrix.
    Alignment: 3xM char matrix showing the alignment results of SeqA and
      SeqB. '|' in the 2nd char row marks locations of matches, while ' '
      in the 1st and 3rd row are unmatched letters.

--------------------------------------------------------------------------
  EXAMPLES
    Case1) "TrimSide" edge cleaning
      SeqA = 'TAATAATTAAT'
      SeqB = 'TCCTAATTGGT'
      [Score, StartAt, MatchAt, Alignment] = alignSeqMEX(SeqA, SeqB)
      Score =
            7
           19
      Alignment =
         TAATAATTAAT
         |  |||||  |
         TCCTAATTGGT
      StartAt =
            1
            1
      MatchAt =
            1
           11

      [Score, ~, ~, Alignment] = alignSeqMEX(SeqA, SeqB, 0, 'n', 'n', 'l')
       Alignment =
         TAATAATTAAT
            |||||  |
         TCCTAATTGGT

      [Score, ~, ~, Alignment] = alignSeqMEX(SeqA, SeqB, 0, 'n', 'n', 'r')
       Alignment =
         TAATAATTAAT
         |  |||||   
         TCCTAATTGGT

      [Score, ~, ~, Alignment] = alignSeqMEX(SeqA, SeqB, 0, 'n', 'n', 'b')
       Alignment =
         TAATAATTAAT
            |||||   
         TCCTAATTGGT

    Case2) "PreferSide" tie-breaking
      SeqA = 'AAAA'
      SeqB = 'CAAAATTAAAATTAAAAC'
      [~, ~, ~, Alignment] = alignSeqMEX(SeqA, SeqB, 0, 'n', 'n', 'n', 'n', 'l')
       Alignment =
         -AAAA-------------
          ||||             
         CAAAATTAAAATTAAAAC

      [~, ~, ~, Alignment] = alignSeqMEX(SeqA, SeqB, 0, 'n', 'n', 'n', 'n', 'r')
       Alignment =
         -------------AAAA-
                      |||| 
         CAAAATTAAAATTAAAAC

      [~, ~, ~, Alignment] = alignSeqMEX(SeqA, SeqB, 0, 'n', 'n', 'n', 'n', 'n')
       Alignment =
         -------AAAA-------
                ||||       
         CAAAATTAAAATTAAAAC

    Case3) Wildcard matching and exact matching
      SeqA = 'CGAAXCAA'
      SeqB = 'ACGAACGA'
      [~, ~, ~, Alignment] = alignSeqMEX(SeqA, SeqB, 0, 'n', 'n')
       Alignment =
         -CGAAXCAA
          ||||| |  
         ACGAACGA-

      [~, ~, ~, Alignment] = alignSeqMEX(SeqA, SeqB, 0, 'n', 'y')
       Alignment =
         CGAAXCAA
            ||| |
         ACGAACGA

    Case4) AllowedMiss scoring changes
      SeqA = 'ACGTGGTA'
      SeqB = 'ACATGATA'
      [Score, ~, ~, Alignment] = alignSeqMEX(SeqA, SeqB, 0, 'n', 'y')
       Score =
           6
          10
       Alignment =
         ACGTGGTA
         || || ||
         ACATGATA

      [Score, ~, ~, Alignment] = alignSeqMEX(SeqA, SeqB, 1/8, 'n', 'y')
        Score =
            6
           19
 
      [Score, ~, ~, Alignment] = alignSeqMEX(SeqA, SeqB, 1/8, 'n', 'y')
        Score =
            6
           19
*/

#include "AlignTool.hpp"
        
void mexFunction(int nlhs,        mxArray *plhs[],
                 int nrhs, const  mxArray *prhs[]) {
    
    if (nrhs < 2 || nrhs > 8) {
        mexErrMsgIdAndTxt("alignSeqMEX:nrhs", "Incorrect number of inputs. Min is 2. Max is 8.");\
    }
    if (nlhs > 4) {
        mexErrMsgIdAndTxt("alignSeqMEX:nlhs", "Too many outputs. Max is 4.");
    }
    
    double MissRate    = 0;
    mxChar Alphabet    = 'n';
    mxChar ExactMatch  = 'n';
    mxChar TrimSide    = 'n';
    mxChar PenaltySide = 'n';
    mxChar PreferSide  = 'n';
    if (nrhs >= 3) {
        if (!mxIsDouble(prhs[2])) { mexErrMsgIdAndTxt("alignSeqMEX:prhs", "Input3: MissRate must be a scalar between 0.0 to 1.0."); }
        MissRate = mxGetScalar(prhs[2]);
        if (nrhs >= 4) {
            if (!mxIsChar(prhs[3])) { mexErrMsgIdAndTxt("alignSeqMEX:prhs", "Input4: Alphabet must be a char 'n' (nucleotide), 'a' (amino acid), or other (basic char)."); }
            Alphabet = *mxGetChars(prhs[3]);
            if (nrhs >= 5) {
                if(!mxIsChar(prhs[4])) { mexErrMsgIdAndTxt("alignSeqMEX:prhs", "Input5: ExactMatch must be char 'n' or 'y'."); }
                ExactMatch = *mxGetChars(prhs[4]);
                if (nrhs >= 6) {
                    if (!mxIsChar(prhs[5])) { mexErrMsgIdAndTxt("alignSeqMEX:prhs", "Input6: TrimSide must be char 'n' (none), 'l' (left), 'r' (right), 'b' (both)."); }
                    TrimSide = *mxGetChars(prhs[5]);
                    if (nrhs >= 7) {
                        if (!mxIsChar(prhs[6])) { mexErrMsgIdAndTxt("alignSeqMEX:prhs", "Input7: PenaltySide must be a char 'n' (none), 'l' (left), 'r' (right), 'b' (both)."); }
                        PenaltySide = *mxGetChars(prhs[6]);
                        if (nrhs >= 8) {
                            if(!mxIsChar(prhs[7])) { mexErrMsgIdAndTxt("alignSeqMEX:prhs", "Input8: PreferSide must be a char 'n' (none), 'l' (left), 'r' (right)."); }
                            PreferSide = *mxGetChars(prhs[7]);
                        }
                    }
                }
            }
        }
    }
    
    mxChar *pSeqA, *pSeqB;
    mwSize LenA, LenB;
    
    if (!mxIsChar(prhs[0]) || mxGetM(prhs[0]) > 1) {
        mexErrMsgIdAndTxt("alignSeqMEX:prhs", "Input1: SeqA must be a 1xN char array.");
    } else {
        pSeqA = mxGetChars(prhs[0]);
        LenA  = mxGetN(prhs[0]);
    }
    
    mwSize Z = mxIsCell(prhs[1]) ? mxGetNumberOfElements(prhs[1]) : 1;
    align_info AI[Z];
    if (mxIsCell(prhs[1])) {
        mwSize MaxN = 0, CurN = 0;
        for (mwSize j = 0; j < Z; j++) {
            CurN = mxGetN(mxGetCell(prhs[1], j));
            if (CurN > MaxN) { MaxN = CurN; }
        }
        bool pMatch[LenA + MaxN];
        for (mwSize j = 0; j < Z; j++) {
            pSeqB = mxGetChars(mxGetCell(prhs[1], j));
            LenB = mxGetN(mxGetCell(prhs[1], j));
            alignSeq(pSeqA, pSeqB, LenA, LenB, MissRate, Alphabet, ExactMatch, TrimSide, PenaltySide, PreferSide, AI[j], pMatch); 
        }
    } else {
        pSeqB = mxGetChars(prhs[1]);
        LenB = mxGetN(prhs[1]);
        bool pMatch[LenA + LenB];
        alignSeq(pSeqA, pSeqB, LenA, LenB, MissRate, Alphabet, ExactMatch, TrimSide, PenaltySide, PreferSide, AI[0], pMatch); 
    }

    int q = 0; //generic counter
    if (nlhs >= 1) {
        plhs[0] = mxCreateDoubleMatrix(2, Z, mxREAL);
        double *pScore = mxGetPr(plhs[0]);
        for (mwSize z = 0; z < Z; z++) {
            pScore[2*z]   = AI[z].Match;
            pScore[2*z+1] = AI[z].Score;
        }
    }

    if (nlhs >= 2) {
        plhs[1] = mxCreateDoubleMatrix(2, Z, mxREAL);
        double *pStartAt = mxGetPr(plhs[1]);
        for (mwSize z = 0; z < Z; z++) {
            pStartAt[2*z]   = 1;
            pStartAt[2*z+1] = AI[z].BShift >= 0 ? AI[z].BShift + 1 : AI[z].BShift;
        }
    }

    if (nlhs >= 3) {
        plhs[2] = mxCreateDoubleMatrix(2, Z, mxREAL);
        double *pMatchAt = mxGetPr(plhs[2]);
        for (mwSize z = 0; z < Z; z++) {
            pMatchAt[2*z]   = AI[z].MatchS + 1;
            pMatchAt[2*z+1] = AI[z].MatchE + 1;
        }
    }

    if (nlhs >= 4) { //Make the alignment
        if (mxIsCell(prhs[1])) {
            plhs[3] = mxCreateCellArray(2, mxGetDimensions(prhs[1]));
            for (mwSize z = 0; z < Z; z++) {
                pSeqB = mxGetChars(mxGetCell(prhs[1], z));
                LenB = mxGetN(mxGetCell(prhs[1], z));
                mxSetCell(plhs[3], z, buildAlignment(pSeqA, pSeqB, LenA, LenB, AI[z]));
            }
        } else {
            plhs[3] = buildAlignment(pSeqA, pSeqB, LenA, LenB, AI[0]);
        }
    }
}