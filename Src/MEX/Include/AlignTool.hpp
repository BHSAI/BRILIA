#ifndef ALIGN_TOOL_HPP
#define ALIGN_TOOL_HPP

#include "mex.h"
#include <math.h>

// align_info will store information of the sequence alignment
struct align_info {
    double Match = 0; //Match: # of nts matched
    double Score = 0; //Score: alignment score
    int BShift = 0;   //BShift: how many letters to shift SeqB right or left of SeqA for alignnment
    int MatchS = -1;  //MatchS: 1st  position where SeqA-to-SeqB alignment position matches (1-based index)
    int MatchE = -1;  //MatchE: last position where SeqA-to-SeqB alignment position matches (1-based index)
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