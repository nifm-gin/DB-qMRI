// nuft_forw.cpp
//
// Function that performs the k-space measurements of a 2D image
//
//                 FOVx-1  FOVy-1
//       m[n] =     sum     sum   I[px,py]*exp(j*2*pi*(kx[n]*px/FOVx+ky[n]*py/FOVy)),
//                 px=0    py=0
// with 0 <= n <= Ns-1.
//
// Inputs: the 2D image I, the 2D k-space trajectory: kx, ky
//
// Output: the vector of k-space measurements.
//
// No gridding is involved. Note that FFT-based implementations can be much faster.
// The implementation aims at minimizing the amount of trigonometric function calls.
// Observed x20 speedup compared to the trivial MATLAB implementation:
//
//	>>m = zeros(1,length(kx));
//	>>[X,Y] = ndgrid((0:size(I,1)-1)/size(I,1),(0:size(I,2)-1)/size(I,2));
//	>>for i=1:length(kx)
//	>>	m(i) = I(:).'*exp(1i*2*pi*(xj(i)*X(:)+yj(i)*Y(:)));
//	>>end
//
// This is a MEX-file for MATLAB.
//
// Copyright, Matthieu Guerquin-Kern, 2012

//#define FAST // Uncomment this if you want the simpler/more reliable/slower code

#include "mex.h"
#include "nuft.h"

/* Input Arguments */

#define	X_IN       prhs[0]
#define	kx_IN      prhs[1]
#define	ky_IN      prhs[2]

/* Output Arguments */

#define	m_OUT   plhs[0]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *X_real,*X_imag,*kx,*ky;
    double *m_real,*m_imag;
    mwSize m,n,Nb_samples;
    
    /* Check for proper number of arguments */
    if (nrhs != 3) {
        mexErrMsgTxt("Three inputs required.");
    }
    else if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }
    
    /* Get the dimensions of X and k */
    m = mxGetM(X_IN);
    n = mxGetN(X_IN);
    Nb_samples = MAX(mxGetM(kx_IN),mxGetN(kx_IN));
    
    /* Create a matrix for the return argument */
    m_OUT = mxCreateDoubleMatrix(1, Nb_samples, mxCOMPLEX);
    
    /* Assign pointers to the various parameters */
    m_real = mxGetPr(m_OUT);
    m_imag = mxGetPi(m_OUT);
    X_real = mxGetPr(X_IN);
    X_imag = mxGetPi(X_IN); // NULL if the matrix is purely real
    kx     = mxGetPr(kx_IN);
    ky     = mxGetPr(ky_IN);
    
    /* Do the actual computations in a subroutine */
    if (X_imag==NULL){
        H_real(X_real,kx,ky,m_real,m_imag,m,n,Nb_samples);
    }else{
        H_complex(X_real,X_imag,kx,ky,m_real,m_imag,m,n,Nb_samples);
    }
    
    return;
}
