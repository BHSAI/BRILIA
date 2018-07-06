/*
convStr2NumMEX will convert a string to a number or a group of cell with 
string to a group of cell with number/string

  Num = convStr2NumMEX(Str)

  INPUT
    Str: a string of numbers separated by a non-decimal character
         Note: Can accept scientific notations
 
  OUTPUT
    Num: a scalar, vector, or []

  EXAMPLE
    Str = '[123 -123 1.23e2 -1.23e2 1.23e-2]';
    Num = convStr2NumMEX(Str)
    Num =
      123.0000 -123.0000  123.0000 -123.0000    0.0123
    
    Num = convStr2NumMEX('34')
    Num = 
       34
*/


#include "mex.h"
#include <string>
#include <math.h>

bool isNumber(mxChar Letter) {
    switch (Letter) {
        case '0': return true;
        case '1': return true;
        case '2': return true;
        case '3': return true;
        case '4': return true;
        case '5': return true;
        case '6': return true;
        case '7': return true;
        case '8': return true;
        case '9': return true;
        default: return false;
    }
}

bool isThereAnyNumberLeft(mxChar *pStr, mwSize Len, int Start) {
    for (int j = Start ; j < Len; j++) {
        if (isNumber(pStr[j])) {
            return true;
        }
    }
    return false;
}

double convCharToDouble(mxChar Letter) {
    switch (Letter) {
        case '0': return 0;
        case '1': return 1;
        case '2': return 2;
        case '3': return 3;
        case '4': return 4;
        case '5': return 5;
        case '6': return 6;
        case '7': return 7;
        case '8': return 8;
        case '9': return 9;
        default: mexErrMsgIdAndTxt("convStr2NumMEX:input", "convCharToDouble: Error processing char to double");
    }
}

int findStartOfNum(mxChar *pStr, mwSize Len, int Start) {
    for (int j = Start; j < Len; j++) {
        if ( isNumber(pStr[j]) || 
             (pStr[j] == '-' && j+1 < Len && isNumber(pStr[j+1])) ) {
            return j;
        }
    }
    return -1;
}

int findEndOfNum(mxChar *pStr, mwSize Len, int Start) {
    if (Start < 0) { return -1; }
    int j = Start + 1;
    for (j; j < Len; j++) {
        if ( (!isNumber(pStr[j]) && pStr[j] != '.') || //Do not end on a decimal point. (EX: "123.4545" do not stop at decimal)
             (j == Len-1 && pStr[j] == '.') ) {        //If the last char is a decimal, end at the prior char. (EX: "123." end at 2, not 3)
            break;
        }
    }
    return j--;
}

int findDecimal(mxChar *pStr, int Start, int End) {
    if (Start < 0 || End < 0) { return -1; }
    for (int j = Start; j < End; j++) {
        if (pStr[j] == '.') { return j; }
    }
    return End++;
}

double convStrToDouble(mxChar *pStr, mwSize Len, int *pStart) {
    int Start = findStartOfNum(pStr, Len, *pStart);
    int End = findEndOfNum(pStr, Len, Start);
    int Decimal = findDecimal(pStr, Start, End);
    
    //Determine the sign of the number
    double Value = 0;
    double Sign = 1;
    if (pStr[Start] == '-') {
        Sign = -1;
        Start++;
    }
    
    //Do the left of decimal
    double Factor = 1;
    for (int k = Decimal-1; k >= Start; k--) {
        Value += Factor * convCharToDouble(pStr[k]);
        Factor *= 10;
    }

    //Do the right of decimal
    Factor = 0.1;
    for (int k = Decimal+1; k < End; k++) {
        Value += Factor * convCharToDouble(pStr[k]);
        Factor /= 10;
    }

    *pStart = End++; //Increment the start index as the next one.
    Value = Sign * Value;
    //Determine if its a scientific notation
    if (Len > End-1 && (pStr[End-1] == 'E' || pStr[End-1] == 'e')) { //attempt scientific notation
        bool NegExp;
        if (pStr[End] == '-') {
            NegExp = true;
            End++;
        } else {
            NegExp = false;
        }

        if (isNumber(pStr[End])) {
            int NewEnd = findEndOfNum(pStr, Len, End);
            double Factor = 1;
            double ExpVal = 0;
            for (int j = NewEnd-1; j >= End; j--) {
                ExpVal += Factor * convCharToDouble(pStr[j]);
                Factor *= 10;
            }
            if (NegExp) {
                Value = Value / pow(10, ExpVal);
            } else {
                Value = Value * pow(10, ExpVal);
            }
            *pStart = NewEnd++; //Adjust for the new end
        }
    }
    
    return Value;
}

mxArray *convStrToMatrix(mxChar *pStr, mwSize Len) {
    double Values[Len];
    int k = 0;
    int Start = 0;
    while (isThereAnyNumberLeft(pStr, Len, Start)) {
        Values[k++] = convStrToDouble(pStr, Len, &Start);
    }
    mxArray *pmxArray = mxCreateDoubleMatrix(1, k, mxREAL);
    double *pOut = mxGetPr(pmxArray); 
    for (mwSize j = 0; j < k; j++) {
        pOut[j] = Values[j];
    }
    return pmxArray;
}

void mexFunction(int nlhs,        mxArray *plhs[],
                 int nrhs, const  mxArray *prhs[]) {
    
    if (nrhs != 1) {
        mexErrMsgIdAndTxt("convStr2NumMEX:input", "Can accept only 1 input.");
    }
    
    if (nrhs == 1 && !(mxIsChar(prhs[0]) || mxIsCell(prhs[0]))) {
        mexErrMsgIdAndTxt("convStr2NumMEX:input", "Input must be a char or cell array.");
    }
    
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("convStr2NumMEX:output", "Too many outputs - max is 1.");
    }
        
    if (mxIsCell(prhs[0])) {
        mwSize M = mxGetM(prhs[0]);
        mwSize N = mxGetN(prhs[0]);
        plhs[0] = mxCreateCellMatrix(M, N);
        for (mwSize j = 0; j < M*N; j++) {
            mwSize Len = mxGetN(mxGetCell(prhs[0], j));
            if (!mxIsChar(mxGetCell(prhs[0], j))) { continue; }
            mxChar *pChar = mxGetChars(mxGetCell(prhs[0], j));
            mxSetCell(plhs[0], j, convStrToMatrix(pChar, Len));
        }
    } else {
        mwSize Len = mxGetN(prhs[0]);
        mxChar *pChar = mxGetChars(prhs[0]);
        plhs[0] = convStrToMatrix(pChar, Len);
    }
}