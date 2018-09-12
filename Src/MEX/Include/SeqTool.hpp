#ifndef SEQ_TOOL_HPP
#define SEQ_TOOL_HPP

#include "mex.h"
#include <string>

int nt2int(mxChar); 
mxChar int2nt(int);
void fixSeqDNA(char*);
void fixSeqRNA(char*);

#endif