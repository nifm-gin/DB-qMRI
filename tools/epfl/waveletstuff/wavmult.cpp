// wavmult.cpp
//
// Function that multiplies a wavelet object with a real constant or a vector of real constants
//
// Matthieu Guerquin-Kern, Biomedical Imaging Group - EPF Lausanne, 2009-16-11

//#include <stdlib.h>
//#include <math.h>
#include "mex.h"
//#include "matrix.h"

/* Input Arguments */

#define	alpha_IN      prhs[0]
#define	w_IN      prhs[1]

/* Output Arguments */

#define	w_OUT   plhs[0]

void multiply(double alpha, double *a, double *c, int Nelmts) {
    for(int i=0;i<Nelmts;i++){
        c[i] = alpha*a[i];
    }
    return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    /* Check for proper number of arguments */
    
    if (nrhs < 2) {
        mexErrMsgTxt("Almost two inputs required.");
    }
    else if (nrhs >= 3){
        mexErrMsgTxt("Too many inputs.");
    }
    
    if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }
    
    int J = mxGetNumberOfElements(w_IN);
    mwSize ndim;
    const mwSize * dims;
    
    int Nalpha = mxGetNumberOfElements(alpha_IN);
    double *alpha = mxGetPr(alpha_IN);
    
    if (w_OUT==NULL){
        //mexPrintf("Allocate space for output\n");
        ndim = mxGetNumberOfDimensions(w_IN);
        dims = mxGetDimensions(w_IN);
        w_OUT = mxCreateCellArray(ndim, dims);
    }else{
        if (J!=mxGetNumberOfElements(w_OUT)){
            mexErrMsgTxt("Output object is of different size.");
        }
    }
    
    //mexPrintf("J: %d\n", J);
    int Nsub, Nelmt, counter=0;
    mxArray *suba, *subc;
    mxArray *elmta, *elmtc;
    double* tabc_r, *tabc_i;
    for(int i=0;i<J;i++){
        //mexPrintf("level: %d\n", i+1);
        suba = mxGetCell(w_IN, i);
        subc = mxGetCell(w_OUT, i);
        Nsub = mxGetNumberOfElements(suba);
        if (subc==NULL){
            //mexPrintf("Allocate space for output\n");
            ndim = mxGetNumberOfDimensions(suba);
            dims = mxGetDimensions(suba);
            subc = mxCreateCellArray(ndim, dims);
        }
        for(int s=0;s<Nsub;s++){
            //mexPrintf("subband: %d\n", s+1);
            elmta = mxGetCell(suba, s);
            elmtc = mxGetCell(subc, s);
            if ((elmta==NULL)||(elmta==NULL)){
                mexErrMsgTxt("Error.");
            }
            Nelmt = mxGetNumberOfElements(elmta);
            if (elmtc==NULL){
                //mexPrintf("Allocate space for output\n");
                ndim = mxGetNumberOfDimensions(elmta);
                dims = mxGetDimensions(elmta);
                if (mxIsComplex(elmta)){
                    elmtc = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxCOMPLEX);
                }else{
                    elmtc = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
                }
            }
            
            if (mxIsComplex(elmta)){
                tabc_i = mxGetPi(elmtc);
                if (Nalpha==1){multiply(alpha[0], mxGetPi(elmta), tabc_i, Nelmt);}
                else if(Nalpha==J+1){
                    multiply(alpha[i], mxGetPi(elmta), tabc_i, Nelmt);
                    //mexPrintf("mult. coeff.: %f\n", alpha[i]);
                }else{
                    multiply(alpha[counter], mxGetPi(elmta), tabc_i, Nelmt);
                    //mexPrintf("counter = %d, mult. coeff.: %f\n", counter,alpha[counter]);
                }
                mxSetPi(elmtc, tabc_i);
            }
            tabc_r = mxGetPr(elmtc);      
            if (Nalpha==1){multiply(alpha[0], mxGetPr(elmta), tabc_r, Nelmt);}
            else if(Nalpha==J+1){
                multiply(alpha[i], mxGetPr(elmta), tabc_r, Nelmt);
                //mexPrintf("mult. coeff.: %f\n", alpha[i]);
            }else{
                multiply(alpha[counter], mxGetPr(elmta), tabc_r, Nelmt);
                //mexPrintf("counter = %d, mult. coeff.: %f\n", counter,alpha[counter]);
            }
            counter++;
            mxSetPr(elmtc, tabc_r);
            mxSetCell(subc , s, elmtc);
        }
        mxSetCell(w_OUT , i, subc);
    }
    return;
}