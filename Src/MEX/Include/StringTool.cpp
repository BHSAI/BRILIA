/*  containsChar will see if a letter matches to anything from the pattern 
 *  char, returning a true or false. This is used mainly used for sequence 
 *  alignment to check if A is a subset of W (A or T).
 *
 *  Result = containsChar(mxChar, std::string)
 *
 *  EXAMPLE
 *    bool Result = containsChar('C', "AT") 
 *   //returns true because 'C' contains one of the char 'A' or 'T'
 */

#include "StringTool.hpp"

bool containsChar(mxChar Letter, std::string Pattern) {
    for (mwSize j = 0; j < Pattern.length(); j++) {
        if (Letter == Pattern[j])
            return true;
    }
    return false;
}