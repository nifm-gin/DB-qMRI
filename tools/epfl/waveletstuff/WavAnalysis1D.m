function xs_hat = WavAnalysis1D(x_hat, htldstr_hat, d, decimation)
% xs_hat = WavAnalysis1D(x_hat, htldstr_hat, d, decimation)
%
% Frequency-domain implementation of the analysis side of a
% single channel within a (dyadically decimated or undecimated)
% filter bank, along the d-th dimension of a signal.
%
% x_hat:      DFT of the signal to be decomposed.
% htldstr:    Complex conjugate of the DFTs of the analysis filter.
% d:          The dimension to be decomposed.
% decimation: True or false.
%
% xs_hat:     DFT of the transform coefficients.
%
% (c) Cedric Vonesch, 2007.10.18-2008.04.04

xs_hat = FilterSep(x_hat, d, htldstr_hat);
if decimation
	xs_hat = DownSampleDim(xs_hat, d);
end