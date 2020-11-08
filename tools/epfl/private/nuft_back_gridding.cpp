// Backward fast Gaussian gridding.
//
// Implements the algorithm by L. Greengard, 2004.
//
// The code handles arbitrary dimensions.
//
// Copyright, Matthieu Guerquin-Kern, 2012

#include "mex.h"
#include "nuft.h"

/* Input Arguments */
#define	V0	prhs[0]
#define	PX  prhs[1]
#define PY  prhs[2]
#define MR  prhs[3]
#define E2X prhs[4]
#define E2Y prhs[5]
#define E3X prhs[6]
#define E3Y prhs[7]

/* Output Arguments */
#define	C	plhs[0]

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    /* Check input and output arguments */
    
    if (nrhs != 8)
        mexErrMsgTxt("Eight input arguments required.");
    if (nlhs > 1)
        mexErrMsgTxt("Too many output arguments.");
    
    if (mxGetNumberOfDimensions(V0)>2)
        mexErrMsgTxt("V0 must have 1 dimension.");
    if (mxGetN(V0)>1&&mxGetM(V0)>1)
        mexErrMsgTxt("V0 must have 1 dimension.");
    if (!mxIsClass(PX, "int32"))
        mexErrMsgTxt("PX is not INT32.");
    if (!mxIsClass(PY, "int32"))
        mexErrMsgTxt("PY is not INT32.");
    if (mxGetN(MR)*mxGetM(MR)!=2)
        mexErrMsgTxt("MR must have 2 elements.");
    /* Assign pointers to the various parameters */
    
    int nb = mxGetN(E2X)*mxGetM(E2X);
    if (mxGetN(E2Y)*mxGetM(E2Y)!=nb)
        mexErrMsgTxt("E2X and E2Y must have the same number of elements.");
    double *mr = (double*)mxGetData(MR);
    int mrx = (int)mr[0];
    int mry = (int)mr[1];
    int msp = mxGetN(E3X)*mxGetM(E3X)-1;
    if (mxGetN(E3Y)*mxGetM(E3Y)!=msp+1)
        mexErrMsgTxt("E3X and E3Y must have the same number of elements.");
    
    double *e2x = mxGetPr(E2X);
    double *e2y = mxGetPr(E2Y);
    double *e3x = mxGetPr(E3X);
    double *e3y = mxGetPr(E3Y);
    
    int *px = (int*)mxGetData(PX);
    if (mxGetN(PX)*mxGetM(PX)!=nb)
        mexErrMsgTxt("PX and E2X must have the same number of elements.");
    int *py = (int*)mxGetData(PY);
    if (mxGetN(PY)*mxGetM(PY)!=nb)
        mexErrMsgTxt("PY and E2X must have the same number of elements.");
    
    if (!mxIsComplex(V0)){
        mexPrintf("\nReal valued input\n");
        C = mxCreateDoubleMatrix(mrx, mry, mxREAL);
        double *v0r = mxGetPr(V0);
        double *cr = mxGetPr(C);
        backward_gaussian_gridding_real(cr,v0r,msp,px,py,mrx,mry,nb,e2x,e2y,e3x,e3y);        
    }else{
        //mexPrintf("\nComplex valued input\n");
        C = mxCreateDoubleMatrix(mrx, mry, mxCOMPLEX);
        double *cr = mxGetPr(C);
        double *ci = mxGetPi(C);
        double *v0r = mxGetPr(V0);
        double *v0i = mxGetPi(V0);
        backward_gaussian_gridding(cr,ci,v0r,v0i,msp,px,py,mrx,mry,nb,e2x,e2y,e3x,e3y);
    }
    
    return;
}
