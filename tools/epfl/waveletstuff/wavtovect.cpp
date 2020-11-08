// wavtovect.cpp
//
// Function that adds two wavelet objects
//
// Matthieu Guerquin-Kern, Biomedical Imaging Group - EPF Lausanne, 2009-13-11

#include <stdlib.h>
//#include <math.h>
#include "mex.h"
//#include "matrix.h"
#include "string.h"

/* Input Arguments */

#define	c_IN      prhs[0]
#define	w_IN      prhs[1]
#define	numel_IN  prhs[2]

/* Output Arguments */

#define	vct_OUT   plhs[0]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* Check for proper number of arguments */
    
    if (nrhs < 3) {
        mexErrMsgTxt("Almost three inputs required.");
    }
    else if (nrhs >= 4){
        mexErrMsgTxt("Too many inputs.");
    }
    
    if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }
    
    int J = mxGetNumberOfElements(w_IN);
    mwSize ndim;
    const mwSize * dims;
    int numel = mxGetPr(numel_IN)[0];
    //mexPrintf("number of elements: %d\n", numel);
    
    vct_OUT = mxCreateDoubleMatrix(numel, 1, mxCOMPLEX);
    
    double *vct_r=mxGetPr(vct_OUT);
    double *vct_i=mxGetPi(vct_OUT);
    
    //mexPrintf("J: %d\n", J);
    int Nsub, shift = 0;
    mwSize sizebuf;
    mxArray *sub, *elmt;
    double *pr = mxGetPr(c_IN);
    double *pi = mxGetPi(c_IN);
    int Nelmt = mxGetNumberOfElements(c_IN);
    if (pi==NULL){
        for(int e=0;e<Nelmt;e++){
            vct_r[shift+e]=pr[e];
            vct_i[shift+e]=(double)0;
        }
    }else{
        for(int e=0;e<Nelmt;e++){
            vct_r[shift+e]=pr[e];
            vct_i[shift+e]=pi[e];
        }
    }
    shift += Nelmt;
    for(int i=0;i<J;i++){
        //mexPrintf("level: %d\n", i+1);
        sub = mxGetCell(w_IN, i);
        Nsub = mxGetNumberOfElements(sub);
        //mexPrintf("Nsub: %d\n", Nsub);
        for(int s=0;s<Nsub;s++){
            elmt = mxGetCell(sub, s);
            sizebuf = mxGetElementSize(elmt);
            Nelmt = mxGetNumberOfElements(elmt);
            //mexPrintf("subband: %d; nb elements: %d\n", s+1,Nelmt);
            pr = mxGetPr(elmt);
            pi = mxGetPi(elmt);
            if (pi==NULL){
                for(int e=0;e<Nelmt;e++){
                    vct_r[shift+e]=pr[e];
                    vct_i[shift+e]=(double)0;
                }
            }else{
                for(int e=0;e<Nelmt;e++){
                    vct_r[shift+e]=pr[e];
                    vct_i[shift+e]=pi[e];
                }
            }
            shift += Nelmt;
            //mexPrintf("Nelements: %d, sizebuf: %d, shift: %d, example: %e vs %e\n", Nelmt, sizebuf, shift, *(vct_r+shift), *pr);
        }
    }
    mxSetPr(vct_OUT, vct_r);
    mxSetPi(vct_OUT, vct_i);
    return;
}