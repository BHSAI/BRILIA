/*  
findPerimeterMEX will get the row, column, and linear index of pixels 
around a certain pixel. This is faster than using MATLAB if/else statements
and sub2ind functions.

  [R, C, Idx] = findPerimeterMEX(Size, Rc, Cc)

  INPUT
    Size: size of the 2D matrix
    Rc: row number of the pixel
    Cc: col number of the pixel   
 
  OUTPUT
    R: row index of the perimeter pixels
    C: column index of the perimeter pixels
    Idx: linear index of the perimeter pixels

  EXAMPLE
    G = rand(10);
    [R, C, Idx] = findPerimeterMEX(size(G), 1, 1)
    R = 
        2
        1
        2
    C = 
        1
        2
        2
    Idx =
        2
       11
       12

*/

#include "mex.h"

void mexFunction(int nlhs,        mxArray *plhs[],
                 int nrhs, const  mxArray *prhs[]) {

    if (nrhs != 3) {
        mexErrMsgIdAndTxt("findPerimeterMEX:nrhs", "Incorrect number of inputs. Expected 3.");
    }
    if (nlhs < 2 || nlhs > 3) {
        mexErrMsgIdAndTxt("findPerimeterMEX:nlhs", "Incorrect number of outputs. Min is 2. Max is 3.");
    }
    if (mxGetM(prhs[0]) != 1 || mxGetN(prhs[0]) != 2) {
        mexErrMsgIdAndTxt("findPerimeterMEX:prhs", "Input must be a 1x2 size matrix.");
    }
    
    
    double *pSize = mxGetPr(prhs[0]);
    double Rc = mxGetScalar(prhs[1]);
    double Cc = mxGetScalar(prhs[2]);

    double R[8];
    double C[8];
    double NewR, NewC;
    int Count = 0;
    for (int c = -1; c <= 1; c++) {
        NewC = Cc + c; 
        if ((NewC < 1) || (NewC > pSize[1])) { continue; }
        for (int r = -1; r <= 1; r++) {
            if (r == 0 && c == 0) { continue; }
            NewR = Rc + r;
            if ((NewR < 1) || (NewR > pSize[0])) { continue; }
            R[Count] = NewR;
            C[Count] = NewC;
            Count++;
        }
    }
    
    plhs[0] = mxCreateDoubleMatrix(Count, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(Count, 1, mxREAL);
    double *pR = mxGetPr(plhs[0]);
    double *pC = mxGetPr(plhs[1]);
    for (int j = 0; j < Count; j++) {
       pR[j] = R[j];
       pC[j] = C[j];
    }

    if (nlhs == 3) {
        plhs[2] = mxCreateDoubleMatrix(Count, 1, mxREAL);
        double *pIdx = mxGetPr(plhs[2]);
        for (int j = 0; j < Count; j++) {
            pIdx[j] = pSize[0]*(pC[j]-1) + pR[j];
        }
    }
}