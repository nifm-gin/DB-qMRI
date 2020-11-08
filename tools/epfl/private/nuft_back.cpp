// nuft_back.cpp
//
// Function that performs the adjoint operation of the k-space measurement of a 2D image
//
//                    Ns
//       I[px,py] =   sum  m[n]*exp(j*2*pi*(kx[n]*px/FOVx+ky[n]*py/FOVy)),
//                    n=1
// with 0 <= px <= FOVx-1, 0 <= py <= FOVy-1.
//
// Inputs: the k-space measurements m, the 2D k-space trajectory: kx,ky ,
// 		and the size of the output image FOVx,FOVy.
//
// Output: a 2D image
//
// This is a MEX-file for MATLAB.
//
// No gridding is involved. Note that FFT-based implementations can be much faster.
// The implementation aims at minimizing the amount of trigonometric function calls.
// Observed x20 speedup compared to the trivial MATLAB implementation:
//
//	>>I = zeros(s);
//	>>[X,Y] = ndgrid((0:s(1)-1)/s(1),(0:s(2)-1)/s(2));
//	>>for i=1:length(kx)
//	>>	I = I+m(i).'*exp(-1i*2*pi*(kx(i)*X+ky(i)*Y));
//	>>end
//
// Copyright, Matthieu Guerquin-Kern, 2012

//#define FAST // Uncomment this if you want the simpler/more reliable/slower code

#include "mex.h"
#include "nuft.h"

/* Input Arguments */

#define	m_IN       prhs[0]
#define	kx_IN      prhs[1]
#define	ky_IN      prhs[2]
#define	size_IN    prhs[3]

/* Output Arguments */

#define	X_OUT   plhs[0]


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *X_real,*X_imag,*kx,*ky;
    double *m_real,*m_imag;
    mwSize m,n,Nb_samples;
    double *size;
    
    /* Check for proper number of arguments */
    
    if (nrhs != 4) {
        mexErrMsgTxt("Four inputs required.");
    }
    else if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }
    
    /* Get the dimensions of X and k */
    /*if (!mxIsInt32(size_IN)) {
     * mexErrMsgTxt("The size vector should contain unsigned 8-bits integers!");
     * }*/
    
    size = mxGetPr(size_IN);
    m = (mwSize)size[0];
    n = (mwSize)size[1];
    Nb_samples = MAX(mxGetM(kx_IN),mxGetN(kx_IN));
    
    //mexPrintf("m:  %d; n: %d; nb_samp: %d\n",m,n,Nb_samples);
    
    /* Create a matrix for the return argument */
    X_OUT = mxCreateDoubleMatrix(m, n, mxCOMPLEX);
    
    /* Assign pointers to the various parameters */
    
    m_real = mxGetPr(m_IN);
    m_imag = mxGetPi(m_IN);
    
    X_real = mxGetPr(X_OUT);
    X_imag = mxGetPi(X_OUT);
    
    kx     = mxGetPr(kx_IN);
    ky     = mxGetPr(ky_IN);
    
    /* Do the actual computations in a subroutine */
    HH(m_real,m_imag,kx,ky,X_real,X_imag,m,n,Nb_samples);
    
    return;
}
