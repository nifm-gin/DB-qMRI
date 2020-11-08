// y_hat = FilterSep(x_hat, df, f_hat)
//
// Matlab MEX-file for fast multidimensional separable filtering.
//
// (c) Cedric Vonesch, 2006.11.06-2006.11.09

//#include <stdlib.h>
#include <math.h>
#include "mex.h"
//#include "matrix.h"

// Global variables
int nxmax, nfmax;
int nslicemax, nsliceend;
	
__inline void multiplyRR(double *x_hat, double *f_hat, double *y_hat) {
	int nx, nf;

	nx = 0;
	if (nslicemax == 1) { // We are filtering dimension 0 (df = 0)
		do {
			nf = 0;
			do {
				y_hat[nx] = f_hat[nf] * x_hat[nx];
				nf++;
				nx++;
			} while (nf<nfmax);
		} while (nx<nxmax);
	} else { // df > 0
		int nsliceend;

		do {
			nf = 0;
			do {
				nsliceend = nx + nslicemax;
				do {
					y_hat[nx] = x_hat[nx] * f_hat[nf];
					nx++;
				} while (nx<nsliceend);
				nf++;
			} while (nf<nfmax);
		} while (nx<nxmax);
	}
}

__inline void multiplyRC(double *x_hat, double *f_hatr, double *f_hati, double *y_hatr, double *y_hati) {
	int nx, nf;

	nx = 0;
	if (nslicemax == 1) { // We are filtering dimension 0 (df = 0)
		do {
			nf = 0;
			do {
				y_hatr[nx] = x_hat[nx] * f_hatr[nf];
				y_hati[nx] = x_hat[nx] * f_hati[nf];
				nf++;
				nx++;
			} while (nf<nfmax);
		} while (nx<nxmax);
	} else { // df > 0
		int nsliceend;

		do {
			nf = 0;
			do {
				nsliceend = nx + nslicemax;
				do {
					y_hatr[nx] = x_hat[nx] * f_hatr[nf];
					y_hati[nx] = x_hat[nx] * f_hati[nf];
					nx++;
				} while (nx<nsliceend);
				nf++;
			} while (nf<nfmax);
		} while (nx<nxmax);
	}
}

__inline void multiplyCR(double *x_hatr, double *x_hati, double *f_hat, double *y_hatr, double *y_hati) {
	int nx, nf;

	nx = 0;
	if (nslicemax == 1) { // We are filtering dimension 0 (df = 0)
		do {
			nf = 0;
			do {
				y_hatr[nx] = x_hatr[nx] * f_hat[nf];
				y_hati[nx] = x_hati[nx] * f_hat[nf];
				nf++;
				nx++;
			} while (nf<nfmax);
		} while (nx<nxmax);
	} else { // df > 0
		int nsliceend;

		do {
			nf = 0;
			do {
				nsliceend = nx + nslicemax;
				do {
					y_hatr[nx] = x_hatr[nx] * f_hat[nf];
					y_hati[nx] = x_hati[nx] * f_hat[nf];
					nx++;
				} while (nx<nsliceend);
				nf++;
			} while (nf<nfmax);
		} while (nx<nxmax);
	}
}

__inline void multiplyCC(double *x_hatr, double *x_hati, double *f_hatr, double *f_hati, double *y_hatr, double *y_hati) {
	int nx, nf;

	nx = 0;
	if (nslicemax == 1) { // We are filtering dimension 0 (df = 0)
		do {
			nf = 0;
			do {
				y_hatr[nx] = x_hatr[nx] * f_hatr[nf] - x_hati[nx] * f_hati[nf];
				y_hati[nx] = x_hatr[nx] * f_hati[nf] + x_hati[nx] * f_hatr[nf];
				nf++;
				nx++;
			} while (nf<nfmax);
		} while (nx<nxmax);
	} else { // df > 0
		int nsliceend;

		do {
			nf = 0;
			do {
				nsliceend = nx + nslicemax;
				do {
					y_hatr[nx] = x_hatr[nx] * f_hatr[nf] - x_hati[nx] * f_hati[nf];
					y_hati[nx] = x_hatr[nx] * f_hati[nf] + x_hati[nx] * f_hatr[nf];
					nx++;
				} while (nx<nsliceend);
				nf++;
			} while (nf<nfmax);
		} while (nx<nxmax);
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int dmax, df, d;
	const mwSize *Nx, *Nf;

	// Check for proper number of input arguments
	if (nrhs!=3)  mexErrMsgTxt("There must be exactly 3 input arguments.");
	if (nlhs>1) mexErrMsgTxt("There must be no more than 1 output argument.");

	// Check for proper type of input arguments
	if (!mxIsDouble(prhs[0])) mexErrMsgTxt("x_hat must be a double-precision array.");
	dmax = (int) mxGetNumberOfDimensions(prhs[0]);
	Nx = mxGetDimensions(prhs[0]);
	nxmax = mxGetNumberOfElements(prhs[0]);
	if (nxmax == 0) mexErrMsgTxt("x_hat cannot be empty.");

	if (mxIsComplex(prhs[1]) | (!mxIsDouble(prhs[1])) | (mxGetNumberOfElements(prhs[1]) != 1)) mexErrMsgTxt("df must be a real scalar.");
	df = ((int) mxGetScalar(prhs[1])) - 1;
	if (df < 0 | dmax <= df) mexErrMsgTxt("df must be one of the dimensions of x_hat.");

	if (!mxIsDouble(prhs[2])) mexErrMsgTxt("f_hat must be a double-precision array.");
	Nf = mxGetDimensions(prhs[2]);
	if (mxGetNumberOfDimensions(prhs[2]) > 2 | (Nf[0] != 1 & Nf[1] != 1)) mexErrMsgTxt("f_hat must be a 1D array.");
	nfmax = mxGetNumberOfElements(prhs[2]);
// 	printf("%d\n", nfmax);
// 	printf("%d\n", Nx[df]);
	if (nfmax == 0) mexErrMsgTxt("f_hat cannot be empty.");
	if (nfmax != Nx[df]) mexErrMsgTxt("The number of elements of f_hat must correspond to the dimension of x_hat to be filtered.");

	// Compute "slice length"
	nslicemax = 1; // If df == 0 this is the correct value
	for (d=0; d<df; d++) nslicemax *= Nx[d];

	// Allocate output variable
	if (mxIsComplex(prhs[0]) | mxIsComplex(prhs[2])) plhs[0] = mxCreateNumericArray(dmax, Nx, mxDOUBLE_CLASS, mxCOMPLEX);
	else plhs[0] = mxCreateNumericArray(dmax, Nx, mxDOUBLE_CLASS, mxREAL);

	// Filter x_hat with f_hat
	if (plhs[0]) {
		if (!mxIsComplex(prhs[0]) & !mxIsComplex(prhs[2])) {
			multiplyRR(mxGetPr(prhs[0]), mxGetPr(prhs[2]), mxGetPr(plhs[0]));
		} else if (!mxIsComplex(prhs[0]) & mxIsComplex(prhs[2])) {
			multiplyRC(mxGetPr(prhs[0]), mxGetPr(prhs[2]), mxGetPi(prhs[2]), mxGetPr(plhs[0]), mxGetPi(plhs[0]));
		} else if (mxIsComplex(prhs[0]) & !mxIsComplex(prhs[2])) {
			multiplyCR(mxGetPr(prhs[0]), mxGetPi(prhs[0]), mxGetPr(prhs[2]), mxGetPr(plhs[0]), mxGetPi(plhs[0]));
		} else {
			multiplyCC(mxGetPr(prhs[0]), mxGetPi(prhs[0]), mxGetPr(prhs[2]), mxGetPi(prhs[2]), mxGetPr(plhs[0]), mxGetPi(plhs[0]));
		}
	} else mexErrMsgTxt("Could not allocate memory for output variable.");
}
