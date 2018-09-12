#ifndef ALIGN_TOOL_HPP
#define ALIGN_TOOL_HPP

#include "mex.h"
#include <math.h>

/* align_info will store information of the sequence alignment
 *   Match: # of nts matched
 *   Score: alignment score
 *   BShift: how many letters to shift SeqB right or left of SeqA for alignnment
 *   MatchS: 1st  position where SeqA-to-SeqB alignment position matches (1-based index)
 *   MatchE: last position where SeqA-to-SeqB alignment position matches (1-based index)
 */
struct align_info {
    double Match  = 0;
    double Score  = 0;
    int BShift = 0;
    int MatchS = -1;
    int MatchE = -1;
}; 

int findFirstMatch(bool*, mwSize);
int findLastMatch(bool*, mwSize);
double calcAlignScore(bool*, mwSize, double, mxChar);
void alignSeq(mxChar*, mxChar*, mwSize, mwSize, double, mxChar, mxChar, mxChar, mxChar, mxChar, align_info&);
void alignSeq(mxChar*, mxChar*, mwSize, mwSize, double, mxChar, mxChar, mxChar, mxChar, mxChar, align_info&, bool*);
void cmprSeq(mxChar*, mxChar*, mwSize, mxChar, bool*);
void trimMatchResults(bool*, mwSize, mxChar);
mxArray *buildAlignment(mxChar*, mxChar*, mwSize, mwSize, align_info);
int countMatch(bool*, mwSize);

#endif