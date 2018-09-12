#ifndef HOTSPOT_TOOL_HPP
#define HOTSPOT_TOOL_HPP

#include "mex.h"
#include "SeqTool.hpp"

extern double pSHM_TENDENCY[5][5];
bool isWRC(mxChar*, mwSize, mwSize);
bool isGYW(mxChar*, mwSize, mwSize);
bool isWA(mxChar*, mwSize, mwSize);
bool isTW(mxChar*, mwSize, mwSize);
bool isHotspot(mxChar*, mwSize, mwSize); 
double labelHotspot(mxChar*, mwSize, mwSize);
double countHotspots(mxChar*, mwSize);
double countHotspots(mxChar*, mwSize Len, double*);
void calcSeqShmScore(mxChar*, mxChar*, mwSize, double*);

#endif