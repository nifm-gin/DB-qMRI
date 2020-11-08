// y_hat = DownSampleDim(x_hat, dd)
//
// Matlab MEX-file for downsampling by a factor
// of 2 along a specified dimension.
// Implementation is done in the frequency domain.
//
// (c) Cedric Vonesch, 2006.11.14

// Test data: a = [1 0 0 0; 0 1 0 0; 0 0 0 2; 0 0 2 0]

//#include <stdlib.h>
#include <math.h>
#include "mex.h"
//#include "matrix.h"

// Global variables
int nxmax, nymax;      // Number of elements in the input and output arrays
int nslicemax;         // "Slice length" (number of pixels to be copied contiguously)

// Does the actual copy-paste job
__inline void alias(double *x_hat, double *y_hat) {
	int nx, ny;        // Linear indices for the input and output data
	int nsliceend;     // Index of the last element of the current slice, plus 1

	ny = nx = 0;
	do {
		nsliceend = ny + nslicemax;
		do {
			y_hat[ny] = (x_hat[nx]+x_hat[nx+nslicemax])/2;
			nx++;
			ny++;
		} while (ny<nsliceend);
		nx += nslicemax;
	} while (ny<nymax);
}

// Matlab interface function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int dmax, dd, d;       // Number of dimensions, dimension to downsample, dimension index
	const mwSize *Nx;      // Size of each dimension of the input array
	mwSize *Ny;            // Size of each dimension of the output array
	double *x_hat, *y_hat; // Pointers to input and output data

	// Check for proper number of arguments
	if (nrhs!=2)  mexErrMsgTxt("There must be exactly 2 input arguments.");
	if (nlhs>1) mexErrMsgTxt("There must be no more than 1 output argument.");

	// Check for proper type of arguments
	if (!mxIsDouble(prhs[0])) mexErrMsgTxt("x_hat must be a double-precision array.");
	dmax = (int) mxGetNumberOfDimensions(prhs[0]);
	Nx = mxGetDimensions(prhs[0]);
	nxmax = mxGetNumberOfElements(prhs[0]);
	if (nxmax < 2) mexErrMsgTxt("x_hat must have at least 2 elements.");

	if (mxIsComplex(prhs[1]) | !mxIsDouble(prhs[1]) | mxGetNumberOfElements(prhs[1])!=1) mexErrMsgTxt("dd must be a real scalar.");
	dd = ((int) mxGetScalar(prhs[1])) - 1;
	if (dd < 0 | dmax <= dd) mexErrMsgTxt("dd must be one of the dimensions of x_hat.");
	if (Nx[dd] & 1) mexErrMsgTxt("The dimension to be downsampled must have even length.");

	// Prepare dimension variable
	Ny = (mwSize *) malloc(sizeof(mwSize) * dmax);
	for (d=0; d<dmax; d++) Ny[d] = Nx[d] / (1 + (d == dd));

	// Compute "slice length"
	nslicemax = 1;
	for (d=0; d<=dd; d++) nslicemax *= Ny[d];
	
	// Allocate output variable
	if (mxIsComplex(prhs[0])) plhs[0] = mxCreateNumericArray(dmax, Ny, mxDOUBLE_CLASS, mxCOMPLEX);
	else plhs[0] = mxCreateNumericArray(dmax, Ny, mxDOUBLE_CLASS, mxREAL);
	
	// Copy-paste
	if (plhs[0]) {
		nymax = nxmax / 2;
		x_hat = mxGetPr(prhs[0]);
		y_hat = mxGetPr(plhs[0]);
		alias(x_hat, y_hat);
		if (mxIsComplex(prhs[0])) {
			x_hat = mxGetPi(prhs[0]);
			y_hat = mxGetPi(plhs[0]);
			alias(x_hat, y_hat);
		}
	}

	// Free dimension variable
	free(Ny);
}
