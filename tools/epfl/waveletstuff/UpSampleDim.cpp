// y_hat = UpSampleDim(x_hat, du)
//
// Matlab MEX-file for upsampling by a factor
// of 2 along a specified dimension.
// Implementation is done in the frequency domain.
//
// (c) Cedric Vonesch, 2006.11.14-2006.11.15

// Test data: a = [1 0 0 0; 0 1 0 0; 0 0 0 2; 0 0 2 0]

//#include <stdlib.h>
#include <math.h>
#include "mex.h"
//#include "matrix.h"

// Global variables
int nxmax, nymax;      // Number of elements in the input and output arrays
int nslicemax;         // "Slice length" (number of pixels to be copied contiguously)

// Does the actual copy-paste job (should be faster)
__inline void periodize(double *x_hat, double *y_hat) {
	int nx, ny;        // Linear indices for the input and output data
	int nsliceend;     // Index of the last element of the current slice, plus 1

	ny = nx = 0;
	do {
		nsliceend = nx + nslicemax;
		do {
			y_hat[ny+nslicemax] = y_hat[ny] = x_hat[nx];
			nx++;
			ny++;
		} while (nx<nsliceend);
		ny += nslicemax;
	} while (ny<nymax);
}

// // Does the actual copy-paste job (slower)
// __inline void periodize(double *x_hat, double *y_hat) {
// 	int nx, ny;        // Linear indices for the input and output data
// 	int nsliceend;     // Index of the last element of the current slice, plus 1
// 
// 	ny = nx = 0;
// 	do {
// 		nsliceend = nx + nslicemax;
// 		do {
// 			y_hat[ny++] = x_hat[nx++];
// 		} while (nx<nsliceend);
// 		nx -= nslicemax;
// 		do {
// 			y_hat[ny++] = x_hat[nx++];
// 		} while (nx<nsliceend);
// 	} while (ny<nymax);
// }

// Matlab interface function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int dmax, du, d;       // Number of dimensions, dimension to upsample, dimension index
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
	if (nxmax < 1) mexErrMsgTxt("x_hat must be non-empty.");
	
	if (mxIsComplex(prhs[1]) | !mxIsDouble(prhs[1]) | mxGetNumberOfElements(prhs[1])!=1) mexErrMsgTxt("du must be a real scalar.");
	du = ((int) mxGetScalar(prhs[1])) - 1;
	if (du < 0 | dmax <= du) mexErrMsgTxt("du must be one of the dimensions of x_hat.");
//	if (Nx[du] & 1) mexErrMsgTxt("The dimension to be upsampled must have even length.");

	// Prepare dimension variable
	Ny = (mwSize *) malloc(sizeof(mwSize) * dmax);
	for (d=0; d<dmax; d++) Ny[d] = (1 + (d == du)) * Nx[d];

	// Compute "slice length"
	nslicemax = 1;
	for (d=0; d<=du; d++) nslicemax *= Nx[d];
	
	// Allocate output variable
	if (mxIsComplex(prhs[0])) plhs[0] = mxCreateNumericArray(dmax, Ny, mxDOUBLE_CLASS, mxCOMPLEX);
	else plhs[0] = mxCreateNumericArray(dmax, Ny, mxDOUBLE_CLASS, mxREAL);
	
	// Copy-paste
	if (plhs[0]) {
		nymax = 2 * nxmax;
		x_hat = mxGetPr(prhs[0]);
		y_hat = mxGetPr(plhs[0]);
		periodize(x_hat, y_hat);
		if (mxIsComplex(prhs[0])) {
			x_hat = mxGetPi(prhs[0]);
			y_hat = mxGetPi(plhs[0]);
			periodize(x_hat, y_hat);
		}
	}

	// Free dimension variable
	free(Ny);
}
