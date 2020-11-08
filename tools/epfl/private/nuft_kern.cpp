// nuft_kern.cpp
//
// Function that computes the convolution kernel related to the operation
// of the k-space measurement of a 2D image followed by its adjoint.
//
//                    Ns
//       G[px,py] =   sum  w[n]*exp(-j*2*pi*(kx[n]*px/FOVx+ky[n]*py/FOVy)),
//                    n=1
// with -FOVx+1 <= px <= FOVx-1 and -FOVy+1 <= py <= FOVy-1.
//
// Inputs: the 2D k-space trajectory, the size of FOV and the optional weighting vector.
//
// Output: the 2D kernel image.
//
// This is a MEX-file for MATLAB.
//
// No gridding is involved. Note that FFT-based implementations can be much faster.
// The implementation aims at minimizing the amount of trigonometric function calls.
// Observed x20 speedup compared to the trivial MATLAB implementation:
//
//	>>G = zeros(2*res-1);
//	>>[X,Y] = ndgrid((-res(1)+1:res(1)-1)/res(1),(-res(2)+1:res(2)-1)/res(2));
//	>>for i=1:length(kx)
//	>>	G = G+w(i)*exp(-1i*2*pi*(kx(i)*X+ky(i)*Y));
//	>>end
//
// Copyright, Matthieu Guerquin-Kern, 2012

//#define FAST // Uncomment this if you want the simpler/more reliable/slower code

#include "mex.h"
#include "nuft.h"

/* Input Arguments */

#define	kx_IN      prhs[0]
#define	ky_IN      prhs[1]
#define	size_IN    prhs[2]
#define	weight_IN  prhs[3]

/* Output Arguments */

#define	G_OUT   plhs[0]

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *G_real,*G_imag,*kx,*ky,*weight;
    mwSize m,n,Nb_samples;
    double *size;
    
    
    /* Check for proper number of arguments */
    
    if (nrhs < 3) {
        mexErrMsgTxt("Almost three inputs required.");
    }
    else if (nrhs >= 5){
        mexErrMsgTxt("Too many inputs.");
    }
    
    if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments.");
    }
    
    size = mxGetPr(size_IN);
    m = (mwSize)size[0];
    n = (mwSize)size[1];
    Nb_samples = MAX(mxGetM(kx_IN),mxGetN(kx_IN));
    
    if (nrhs == 3) {
        weight = ((double*)mxCalloc( Nb_samples, sizeof(double)));
        int ind;
        for (ind=0;ind < Nb_samples; ind++){
            weight[ind] = 1;
        }
    }else{
        if (mxIsComplex(weight_IN)) {
            mexWarnMsgTxt("The weighting vector should be real valued. Taking the real part...");
        }
        weight = mxGetPr(weight_IN);
    }
    
    /* Create a matrix for the return argument */
    G_OUT = mxCreateDoubleMatrix(2*m-1, 2*n-1, mxCOMPLEX);
    
    /* Assign pointers to the various parameters */
    
    G_real = mxGetPr(G_OUT);
    G_imag = mxGetPi(G_OUT);
    
    kx     = mxGetPr(kx_IN);
    ky     = mxGetPr(ky_IN);
    
    /* Do the actual computations in a subroutine */
    G(kx,ky,G_real,G_imag,weight,m,n,Nb_samples);
    
    if (nrhs == 3) {
        mxFree(weight);
    }
    
    return;
}
