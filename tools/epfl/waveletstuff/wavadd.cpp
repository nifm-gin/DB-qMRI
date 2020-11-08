// wavadd.cpp
//
// Function that adds two wavelet objects
//
// Matthieu Guerquin-Kern, Biomedical Imaging Group - EPF Lausanne, 2009-13-11

//#include <stdlib.h>
//#include <math.h>
#include "mex.h"
//#include "matrix.h"

/* Input Arguments */

#define	wa_IN      prhs[0]
#define	wb_IN      prhs[1]

/* Output Arguments */

#define	wc_OUT   plhs[0]

void plus(double *a, double *b, double *c, int Nelmts) {
    //mexPrintf("plus %d %d %d\n",Nelmts,(int)(a == NULL),(int)(b == NULL));
    if (a == NULL){
        for(int i=0;i<Nelmts;i++){
            c[i] = b[i];
        }
        return;
    }
    if (b == NULL){
        for(int i=0;i<Nelmts;i++){
            c[i] = a[i];
        }
        return;
    }
    for(int i=0;i<Nelmts;i++){
        c[i] = a[i] + b[i];
    }
    return;
}

void minus(double *a, double *b, double *c, int Nelmts) {
    //mexPrintf("minus %d %d %d\n",Nelmts,(int)(a == NULL),(int)(b == NULL));
    if (a == NULL){
        for(int i=0;i<Nelmts;i++){
            c[i] = -b[i];
        }
        return;
    }
    if (b == NULL){
        for(int i=0;i<Nelmts;i++){
            c[i] = a[i];
        }
        return;
    }
    for(int i=0;i<Nelmts;i++){
        c[i] = a[i] - b[i];
    }
    return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    /* Check for proper number of arguments */
    
    if (nrhs < 2) {
        mexErrMsgTxt("Almost two inputs required.");
    }
    else if (nrhs >= 4){
        mexErrMsgTxt("Too many inputs.");
    }
    int dir;
    if (nrhs == 3){
        double* pm = mxGetPr(prhs[2]);
        dir = (int)pm[0];
        //mexPrintf("%d\n",dir);
    }else{
        dir= (int)1; // Addition by default
    }
    
    if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }
    
    if ((mxIsCell(wa_IN)==0)||(mxIsCell(wb_IN)==0)) {
        mexErrMsgTxt("Input objects are not only wavelets.");
    }
    
    int J = mxGetNumberOfElements(wa_IN);
    mwSize ndim;
    const mwSize * dims;
    
    if (J!=mxGetNumberOfElements(wb_IN)){
        mexErrMsgTxt("Input objects are of different sizes.");
    }
    
    if (wc_OUT==NULL){
        //mexPrintf("Allocate space for output\n");
        ndim = mxGetNumberOfDimensions(wa_IN);
        dims = mxGetDimensions(wa_IN);
        wc_OUT = mxCreateCellArray(ndim, dims);
    }else{
        if (J!=mxGetNumberOfElements(wc_OUT)){
            mexErrMsgTxt("Output object is of different size.");
        }
    }
    
    //mexPrintf("J: %d\n", J);
    int Nsub, Nelmt;
    mxArray *suba, *subb, *subc;
    mxArray *elmta, *elmtb, *elmtc;
    //mxArray *rhs[2];
    //mxArray *lhs[1];
    double* tabc_r, *tabc_i;
    for(int i=0;i<J;i++){
        //mexPrintf("level: %d\n", i+1);
        suba = mxGetCell(wa_IN, i);
        subb = mxGetCell(wb_IN, i);
        subc = mxGetCell(wc_OUT, i);
        if ((suba==NULL)||(subb==NULL)){
            mexErrMsgTxt("Error.");
        }
        Nsub = mxGetNumberOfElements(suba);
        //mexPrintf("Nsub: %d\n", Nsub);
        if (Nsub!=mxGetNumberOfElements(subb)){
            mexErrMsgTxt("Input objects are of different sizes.");
        }
        if (subc==NULL){
            //mexPrintf("Allocate space for output\n");
            ndim = mxGetNumberOfDimensions(suba);
            dims = mxGetDimensions(suba);
            subc = mxCreateCellArray(ndim, dims);
        }
        for(int s=0;s<Nsub;s++){
            //mexPrintf("subband: %d\n", s+1);
            elmta = mxGetCell(suba, s);
            elmtb = mxGetCell(subb, s);
            elmtc = mxGetCell(subc, s);
            if ((elmta==NULL)||(elmta==NULL)){
                mexErrMsgTxt("Error.");
            }
            Nelmt = mxGetNumberOfElements(elmta);
            if (Nelmt!=mxGetNumberOfElements(elmtb)){
                mexErrMsgTxt("Input objects are of different sizes.");
            }
            if (elmtc==NULL){
                //mexPrintf("Allocate space for output\n");
                ndim = mxGetNumberOfDimensions(elmta);
                dims = mxGetDimensions(elmta);
                if (mxIsComplex(elmta)||mxIsComplex(elmtb)){
                    elmtc = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxCOMPLEX);
                }else{
                    elmtc = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
                }
            }
            //rhs[0]=elmta;
            //rhs[1]=elmtb;
            tabc_r = mxGetPr(elmtc);
            if (mxIsComplex(elmta)||mxIsComplex(elmtb)){tabc_i = mxGetPi(elmtc);}
            
            if (dir==1){
                //mexCallMATLAB(1, lhs, 2, rhs, "plus");
                plus(mxGetPr(elmta), mxGetPr(elmtb), tabc_r, Nelmt);
                if (mxIsComplex(elmta)||mxIsComplex(elmtb)){
                    plus(mxGetPi(elmta), mxGetPi(elmtb), tabc_i, Nelmt);
                }
            }else{
                //mexCallMATLAB(1, lhs, 2, rhs, "minus");
                minus(mxGetPr(elmta), mxGetPr(elmtb), tabc_r, Nelmt);
                if (mxIsComplex(elmta)||mxIsComplex(elmtb)){
                    minus(mxGetPi(elmta), mxGetPi(elmtb), tabc_i, Nelmt);
                }
            }
            mxSetPr(elmtc, tabc_r);
            if (mxIsComplex(elmta)){mxSetPi(elmtc, tabc_i);}
            //mxSetCell(subc , s, lhs[0]);
            mxSetCell(subc , s, elmtc);
        }
        mxSetCell(wc_OUT , i, subc);
    }
    return;
}